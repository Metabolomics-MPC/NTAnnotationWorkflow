#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
#######################################################################
library(shiny)
library(plotly)
library(Spectra)
library(stringr)
library(MetaboAnnotation)
library(tidyverse)
library(DT)
#######################################################################
#' Function for generating a plotly headtail plot for one query in mtch object
plotly_headtail <- function(query_spectrum, target_spectrum){
    
    # create data frames for plotting
    top <- data.frame(mz = unlist(mz(query_spectrum)),
                      int = unlist(intensity(query_spectrum)))
    # create layout
    layout <- list(
        title = "",
        xaxis = list(title = "m/z",
                     zeroline=TRUE,
                     range=c(0, max(top$mz)),
                     nticks=8,
                     autorange = TRUE),
        yaxis = list(title = "Signal Intensity [%]",
                     zeroline=TRUE,
                     tickmode='array',
                     tickvals=c(-100, -50, 0, 50, 100),
                     ticktext=c('100','50', '0', '50', '100'))
    )
    
    # create plot
    p <- plot_ly(
        top,
        x =  ~ mz,
        y =  ~ int,
        showlegend = F,
        type = 'bar',
        marker = list(size = 3, color = 'red'),
        hoverinfo = 'none'
    )
    
    p <-
        add_markers(
            p,
            type = "scatter",
            x = top$mz,
            y = top$int,
            hovertemplate = paste('<br>mz:', '%{x}', '<br>int: %{y}<br>'),
            hoverlabel = list(namelength = 0)
        )
    
    if(length(target_spectrum) > 0) {
        
        bottom <- data.frame(mz = unlist(mz(target_spectrum)),
                             int = unlist(intensity(target_spectrum)))
        
        p <-
            add_trace(
                p,
                type = "bar",
                x = bottom$mz,
                y = -bottom$int,
                marker = list(color = 'blue'),
                hoverinfo = 'none'
            )
        
        p <-
            add_markers(
                p,
                x = bottom$mz,
                y = -bottom$int,
                type = 'scatter',
                marker = list(color = 'blue'),
                hovertemplate = paste('<br>mz:', '%{x}', '<br>int: %{y}<br>'),
                hoverlabel = list(namelength = 0)
            )
    }
    
    p <-
        layout(
            p,
            title = layout$title,
            xaxis = layout$xaxis,
            yaxis = layout$yaxis
        )
    
    p <-
        add_annotations(
            p,
            type = 'text',
            x = c(15, 15),
            y = c(100, -100),
            text = c("query", "library"),
            textfont = list(color = c('red', 'blue')),
            showarrow = F
        )
    
    p <- p %>% layout(hovermode = "x", hoverdistance = 1)
    
    # return plot
    p
}
#######################################################################
load_mtch <- function(file){
    #load matched object
    mtch <- readRDS(file)
    print(mtch)
    #get query matches
    mtch_sub <- mtch[whichQuery(mtch)]
    mtch_sub <- pruneTarget(mtch_sub)
    if(length(mtch_sub)==0){
        message("MatchedSpectra object does not contain any matches.")
    }
    return(mtch_sub)
}
#######################################################################
create_choices <- function(mtch_sub){
    #generate choices name for list on side
    l_choices <- as.list(seq(1, length(mtch_sub)))
    names(l_choices) <- paste0(seq(1, length(mtch_sub)),
                               " - MZ",
                               round(query(mtch_sub)$precursorMz,4),
                               "@RT",
                               round(query(mtch_sub)$rtime/60, 2), " min")
    return(l_choices)
}
#######################################################################
create_boolean <- function(mtch_sub){
    # create boolean values
    boolean_values <- list()
    for(i in 1:length(mtch_sub)) {
        x <- mtch_sub[i]
        target_length <- length(x@target[x@matches$target_idx])
        boolean_values[[i]] <- rep(TRUE, length(x@target[x@matches$target_idx]))
    }
    return(boolean_values)
}
#######################################################################
create_table <- function(x){
    # table with matches
    mtch_tbl <- matchedData(x) %>% 
        as_tibble() %>% 
        select("precursorMz",
               "target_precursorMz",
               "rtime",
               "target_rtime",
               "target_name",
               "score",
               "reverse_score",
               "presence_ratio") %>% 
        mutate(score = round(score, 3),
               reverse_score = round(reverse_score, 3),
               presence_ratio = round(presence_ratio, 3)) %>% 
        rename("Query m/z" = "precursorMz",
               "Target m/z" = "target_precursorMz",
               "Query RT" = "rtime",
               "Target RT" = "target_rtime",
               "Name" = "target_name",
               "Forward Score" = "score",
               "Reverse Score" = "reverse_score",
               "Presence Ratio" = "presence_ratio")
}
#######################################################################
verified_only <- function(mtch_sub, boolean_values){
    
    filter_idx <- mtch_sub@matches[which(unlist(boolean_values)),]
    
    mtch_sub_verified <- filterMatches(mtch_sub,
                                       index = as.integer(row.names(filter_idx)))
}
#######################################################################
# set serial processing (for windows)
register(SerialParam())

#file from argument line
file1 <- args[1]  

# change maximum upload file size
options(shiny.maxRequestSize=200*1024^2)
#######################################################################
# Define UI for application
ui <- fluidPage(
  titlePanel(
    # title with logos
    div("ShinyMetaboAnnotation",
        img(src='hmgu.png', width = "180px", align = "right"),
        img(src='ufz.png', width = "180px", align = "right"),
        class = "pull-right")
  ),
  sidebarLayout(
    # define sidebar with list of matches and true match false match buttons
    sidebarPanel(
      # define File upload
      fileInput("file", 'MachedSpectraObject as *rds:', accept=".rds"),
      # define Feature selection
      selectInput('selection', 'Features:', choices=list(), multiple=FALSE, selectize=FALSE, size=25),
      actionButton("b_store", "Save verification", width='157px'),
      downloadButton("d_store", "Select file location")
    ),
    # define main window with text, plotly plot
    mainPanel(
      plotlyOutput("plot"),
      DT::dataTableOutput("dynamic"),
      radioButtons("veri",
                   label = "Hit correct?",
                   choices = list("Yes" = T, "No" = F), 
                   selected = T)
    )
  )
)
#######################################################################
# Define server logic
server <- function(input, output) {
  # define storage places for reactive value that only change if new object is loaded 
  if(!is.na(file1)){
      mtch_sub <- load_mtch(file1)

      q <- reactiveValues(mtch_sub = mtch_sub,
                          l_choices = create_choices(mtch_sub),
                          boolean_values = create_boolean(mtch_sub))
      observe(updateSelectInput(inputId='selection', choices = q$l_choices))
  }else{
      message("No input delivered by run")
      #make placeholder until file upload is performed
      q <- reactiveValues(mtch_sub = MatchedSpectra(),
                          l_choices = 1L,
                          boolean_values = 1L)
  }

  # define storage places for reactive values that change during the shiny application run
  v <- reactiveValues(i = 1)
  id <- reactiveValues(target_idx = 1L)
  specs <- reactiveValues(match = MatchedSpectra(),
                          query = Spectra(),
                          target = Spectra())
  
  # Define what happens for file-upload
  observeEvent(input$file, {
      f <- input$file[1,4]
      q$mtch_sub <- load_mtch(f)
      q$l_choices <- create_choices(q$mtch_sub)
      q$boolean_values <- create_boolean(q$mtch_sub)
      # update selection list
      updateSelectInput(inputId='selection', choices = q$l_choices)
  })
  
  # Define what happens if a selection in list clicked
  observeEvent(input$selection, {
    
    # change index to selection
    v$i <- as.numeric(input$selection)
    
    # filter MatchedSpectra based on index from selection
    x <- q$mtch_sub[v$i]
    mtch_tbl <- create_table(x)

    # create dynamic table
    output$dynamic <- DT::renderDT(DT::datatable(bind_cols(mtch_tbl,
                                                           Hit = q$boolean_values[[v$i]]),
                                             selection = "single"),
                                   server = TRUE)
    
    id$target_idx <- input$dynamic_rows_selected
    output$row <- renderPrint(input$dynamic_rows_selected)
    
    # update reactiveValues containing spectra
    specs$match <- x
    specs$query <- query(x)
    specs$target <- x@target[x@matches$target_idx][input$dynamic_rows_selected]
    
    # render new plotly plot
    observe(output$plot <-  renderPlotly(plotly_headtail(specs$query, specs$target)))

    # change the button to the previous selected verification value
    updateRadioButtons(inputId="veri",
                       choices = list("Yes" = T, "No" = F), 
                       selected = q$boolean_values[v$i][id$target_idx])
  })
  
  # Define what happens if a row in DT is selected
  observeEvent(input$dynamic_rows_selected, {
    
    # change index to selection
    id$target_idx <- input$dynamic_rows_selected
    
    # filter MatchedSpectra based on index from selection
    v$i <- as.numeric(input$selection)
    x <- q$mtch_sub[v$i]

    # update reactiveValues containing spectra
    specs$target <- x@target[x@matches$target_idx][input$dynamic_rows_selected]

    # render new plotly plot
    observe(output$plot <-  renderPlotly(plotly_headtail(specs$query, specs$target)))
    
    # change the button to the previous selected verification value
    updateRadioButtons(inputId="veri",
                       choices = list("Yes" = T, "No" = F),
                       selected = q$boolean_values[v$i][input$dynamic_rows_selected])
    
  })

  # Define what to do if click on radio button
  observeEvent(input$veri, {
    
    # store value in logical vector
      q$boolean_values[[v$i]][id$target_idx] <- as.logical(input$veri)

  })
  
  # Store Result from commandline
  observeEvent(input$b_store, {
    
    if(is.null(input$file)){
        # save depending on the name given from the commandline input
        saveRDS(verified_only(q$mtch_sub, q$boolean_values), paste0(str_replace(file1, ".rds$", "_verified.rds")))
        message(paste0("Stored under ", str_replace(file1, ".rds$", "_verified.rds")))
    }
    
  })
  
  # Store Result for uploaded files
  output$d_store <- downloadHandler(
      filename = function() {
          paste0("verified.rds")
      },
      content = function(file) {
          # save depending on the name given from the commandline input
          saveRDS(verified_only(q$mtch_sub, q$boolean_values), file)
          message("Verified MatchedSpectra object has beeen stored.")
      }
  )
  
}
#######################################################################
# Run the application 
shinyApp(ui = ui, server = server)
