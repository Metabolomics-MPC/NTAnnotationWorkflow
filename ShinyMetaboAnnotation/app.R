#!/usr/bin/env Rscript
#args = commandArgs(trailingOnly=TRUE)
#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
##########################################################################
# Define which object to load => CHANGE!
#filename <- args[1] #"~/git/MetaboliteAnnotationWorkflow/test_output/Annotation_MS2_inhouse/pos_MassbankRecord_Pos_ms2annotation.rds"
filename <- "E:/04_BGC/Project Data/Bio-Chemoinformatics/R/Projects/MetaboliteAnnotationWorkflow/test_output/Annotation_MS2_inhouse/pos_MassbankRecord_Pos_ms2annotation.rds"
#filename <- "E:/04_BGC/Project Data/Bio-Chemoinformatics/R/Projects/MetaboliteAnnotationWorkflow/test_output/Annotation_MS2_external/pos_MSDIal_MSMS-Public-Pos-VS15_ms2annotation.rds"
##########################################################################
library(shiny)
library(plotly)
library(Spectra)
library(stringr)
library(MetaboAnnotation)
library(tidyverse)
library(DT)
#######################################################################
# define plot function
#' Function for generating a plotly headtail plot for one query in mtch object
# plotly_headtail <- function(x, target_idx = 1L){
plotly_headtail <- function(query_spectrum, target_spectrum){
  
  # # isolate query and target spectra
  # query_spectrum <- query(x)
  # target_spectra <- x@target[x@matches$target_idx]
  # 
  # # get corresponding target spectrum
  # target_spectrum <- target_spectra[target_idx]
  
  # create data frames for plotting
  top <- data.frame(mz = unlist(mz(query_spectrum)),
                    int = unlist(intensity(query_spectrum)))
  top <- top[-which(top[,2]==0),]
  
  bottom <- data.frame(mz = unlist(mz(target_spectrum)),
                       int = unlist(intensity(target_spectrum)))
  
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
  p <- plot_ly(top,x=~mz, y=~int, showlegend=F, type='bar',
               marker = list(size = 3, color= 'red'),  hoverinfo='none')
  p <- add_markers(p, type="scatter", x=top$mz, y=top$int, hovertemplate = paste('<br>mz:', '%{x}', '<br>int: %{y}<br>'), hoverlabel = list(namelength=0))
  p <- add_trace(p, type="bar", x=bottom$mz, y=-bottom$int, marker=list(color='blue'), hoverinfo='none')
  p <- add_markers(p, x=bottom$mz, y=-bottom$int, type='scatter', marker=list(color='blue'), hovertemplate = paste('<br>mz:', '%{x}', '<br>int: %{y}<br>'), hoverlabel = list(namelength=0))
  p <- layout(p, title=layout$title, xaxis=layout$xaxis, yaxis=layout$yaxis)
  p <- add_annotations(p, type='text', x=c(15,15), y=c(100,-100), text=c("query", "library"), textfont=list(color=c('red', 'blue')), showarrow=F)
  p <- p %>%  layout(hovermode = "x", hoverdistance=1)
  
  # return plot
  p
}

#######################################################################
#load matched object
mtch <- readRDS(filename)
#get querry matches
mtch_sub <- mtch[whichQuery(mtch)]
#generate choices name for list on side
l_choices <- as.list(seq(1, length(mtch_sub)))
names(l_choices) <- paste0(seq(1, length(mtch_sub)),
                           " - MZ",
                           round(query(mtch_sub)$precursorMz,4),
                           "@RT",
                           round(query(mtch_sub)$rtime, 0))

# create boolean values
boolean_values <- list()
for(i in 1:length(mtch_sub)) {
  x <- mtch_sub[i]
  target_length <- length(x@target[x@matches$target_idx])
  boolean_values[[i]] <- rep(TRUE, length(x@target[x@matches$target_idx]))
}

print(mtch_sub)
print(boolean_values)
print(unlist(boolean_values))

#######################################################################
# Define UI for application
ui <- fluidPage(
  titlePanel(
    # title with logos
    div("ShinyMetaboAnnotation",
        img(src='hmgu.png', width = "250px", align = "right"),
        img(src='ufz.png', width = "250px", align = "right"),
        class = "pull-right")
  ),
  sidebarLayout(
    # define sidebar with list of matches and true match false match buttons
    sidebarPanel(
      selectInput('selection', 'Features:', choices=l_choices , multiple=TRUE, selectize=FALSE, size=25),
      # actionButton("b_back", "back", width='275px'),
      # actionButton("b_next", "next", width='275px'),
      # actionButton("b_close", "load new object", width='275px'),
      actionButton("b_store", "store verification", width='275px')
    ),
    # define main window with text, plotly plot and buttons
    mainPanel(
      DT::dataTableOutput("dynamic"),
      radioButtons("veri",
                   label = "Hit",
                   choices = list("True match" = T, "False match" = F), 
                   selected = T),
      plotlyOutput("plot")
    )
  )
)
#######################################################################
# Define server logic
server <- function(input, output) {
  # define storage places for reactive values that change during the shiny application run
  v <- reactiveValues(i = 1)
  r <- reactiveValues(veri_logical = boolean_values)
  id <- reactiveValues(target_idx = 1L)
  specs <- reactiveValues(match = MatchedSpectra(),
                          query = Spectra(),
                          target = Spectra())
  
  # selection in list
  observeEvent(input$selection, {
    
    # change index to selection
    v$i <- as.numeric(input$selection)
    
    # filter MatchedSpectra based on index from selection
    x <- mtch_sub[v$i]

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
    
    # create dynamic table
    output$dynamic <- DT::renderDT(DT::datatable(bind_cols(mtch_tbl,
                                                           hit = r$veri_logical[[v$i]]),
                                             selection = "single"),
                                   server = TRUE,
                                   options = list(pageLength = 10))
    
    id$target_idx <- input$dynamic_rows_selected
    output$row <- renderPrint(input$dynamic_rows_selected)
    
    # update reactiveValues containing spectra
    specs$match <- x
    specs$query <- query(x)
    
  })
  
  # row selection
  observeEvent(input$dynamic_rows_selected, {
    
    # change index to selection
    # filter MatchedSpectra based on index from selection
    v$i <- as.numeric(input$selection)
    x <- mtch_sub[v$i]
    id$target_idx <- input$dynamic_rows_selected
    
    # update reactiveValues containing spectra
    specs$target <- x@target[x@matches$target_idx][input$dynamic_rows_selected]
    
    # change the button to the previous selected verification value
    updateRadioButtons(inputId="veri",
                       choices = list("True match" = T, "False match" = F), 
                       selected = r$veri_logical[v$i][input$dynamic_rows_selected])
    
    # render new plotly plot
    # observe(output$plot <-  renderPlotly(plotly_headtail(x, id$target_idx)))
    observe(output$plot <-  renderPlotly(plotly_headtail(specs$query, specs$target)))
    
  })

  # Define what to do if click on radio button
  observeEvent(input$veri, {
    
    # store value in logical vector
    r$veri_logical[[v$i]][id$target_idx] <- as.logical(input$veri)

  })
  
  #Store Result
  observeEvent(input$b_store, {
    
    # subset the match object
    print(r$veri_logical)
    filter_idx <- mtch_sub@matches[which(unlist(r$veri_logical)),]
    print(filter_idx)
    print(filterMatches(mtch_sub, index = as.integer(row.names(filter_idx))))
    mtch_sub_verified <- filterMatches(mtch_sub,
                                       index = as.integer(row.names(filter_idx)))
    
    # save as new rds
    saveRDS(mtch_sub_verified, paste0(str_replace(filename, ".rds$", "_verified.rds")))
    
  })
}
#######################################################################
# Run the application 
shinyApp(ui = ui, server = server)
