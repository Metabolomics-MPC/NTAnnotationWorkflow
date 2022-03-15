#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
##########################################################################
# Define which object to load => CHANGE!
filename <- "~/git/MetaboliteAnnotationWorkflow/test_output/Annotation_MS2_inhouse/pos_MassbankRecord_Pos_ms2annotation.rds"
##########################################################################
library(shiny)
library(plotly)
library(Spectra)
library(stringr)
library(MetaboAnnotation)
#######################################################################
# define plot function
#' Function for generating a plotly headtail plot for one query in mtch object
plotly_headtail <- function(x){
    
    top <- data.frame(mz=unlist(mz(query(x))), int=unlist(intensity(query(x))))
    top <- top[-which(top[,2]==0),]
    y <- x@target[x@matches$target_idx]
    
    bottom <- data.frame(mz=unlist(mz(y)),int=unlist(intensity(y)))
    
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
    
    p <- plot_ly(top,x=~mz, y=~int, showlegend=F, type='bar',
                 marker = list(size = 3, color= 'red'),  hoverinfo='none')
    p <- add_markers(p, type="scatter", x=top$mz, y=top$int, hovertemplate = paste('<br>mz:', '%{x}', '<br>int: %{y}<br>'), hoverlabel = list(namelength=0))
    p <- add_trace(p, type="bar", x=bottom$mz, y=-bottom$int, marker=list(color='blue'), hoverinfo='none')
    p <- add_markers(p, x=bottom$mz, y=-bottom$int, type='scatter', marker=list(color='blue'), hovertemplate = paste('<br>mz:', '%{x}', '<br>int: %{y}<br>'), hoverlabel = list(namelength=0))
    p <- layout(p, title=layout$title, xaxis=layout$xaxis, yaxis=layout$yaxis)
    p <- add_annotations(p, type='text', x=c(15,15), y=c(100,-100), text=c("query", "library"), textfont=list(color=c('red', 'blue')), showarrow=F)
    p <- p %>%  layout(hovermode = "x", hoverdistance=1)
    return(p)
}
#######################################################################
#load matched object
mtch <- readRDS(filename)
#get querry matches
mtch_sub <- mtch[whichQuery(mtch)]
#generate choices name for list on side
l_choices <- as.list(seq(1, length(mtch_sub)))
names(l_choices) <- paste(seq(1, length(mtch_sub)),mtch_sub$target_name, "-", rep(T, length(mtch_sub)))
#######################################################################
# Define UI for application
ui <- fluidPage(
    # define sidebar with list of matches and true match false match buttons
    sidebarPanel(
        selectInput('selection', 'matches:', choices=l_choices , multiple=TRUE, selectize=FALSE, size=25),
        radioButtons("veri",
                     label="",
                     choices = list("True match" = T, "False match" = F), 
                     selected=T),
        width=3
    ),
    # define main window with text, plotly plot and buttons
    mainPanel(
    htmlOutput("t_number"),
    htmlOutput("t_formula"),
    htmlOutput("t_query"),
    htmlOutput("t_library"),
    htmlOutput("t_score"),
    plotlyOutput("plot"),
    actionButton("b_back", "back", width='275px'),
    actionButton("b_next", "next", width='275px'),
    actionButton("b_close", "load new object", width='275px'),
    actionButton("b_store", "store verification", width='275px')
    )
)
#######################################################################
# Define server logic
server <- function(input, output) {
    # define storage places for reactive values that change during the shiny application run
    v <- reactiveValues(i = 0)
    r <- reactiveValues(veri_logical=rep(T, length=length(mtch_sub)))
    
    # Go forward with next button
    observeEvent(input$b_next, {
        
        # only react if you are inside the length of the matched object
        if(v$i %in% 0:(length(mtch_sub)-1)){
            
            # increase index
            v$i <- v$i + 1
            
            # chose new match
            x <- mtch_sub[v$i]
            
            # change to new text
            observe(output$t_number <- renderText(HTML(paste0("<b>", v$i, "/", length(mtch_sub), " ", x$target_name, "</b>"))))
            observe(output$t_formula <- renderText(HTML(paste0("Formula: ", x$target_formula, " InchiKey: ", x$target_inchikey))))
            observe(output$t_query <- renderText(HTML(paste0("Query: ", round(x$precursorMz,5), " @ ", round(x$rtime/60,2)))))
            observe(output$t_library <- renderText(HTML(paste0("Library: ", round(x$target_precursorMz,4), " @ ", round(x$target_rtime/60,2)))))
            observe(output$t_score <- renderText(HTML(paste0("Score: ", round(x$score, 3), " Reverse: ", round(x$reverse_score, 3)))))
            
            # render new plotly plot
            output$plot <-  renderPlotly(plotly_headtail(mtch_sub[v$i]))
            
            # change preselection in input table to new value => works, but somehow slows down the process
            #updateSelectInput(session = getDefaultReactiveDomain(), inputId = 'selection', selected=v$i)
            
            # change the button to the previous selected verification value
            updateRadioButtons(inputId="veri",
                               choices = list("True match" = T, "False match" = F), 
                               selected=r$veri_logical[v$i])
        }
    })
    # Go back with button
    observeEvent(input$b_back, {
        # only react if you are inside the length of the matched object
        if(v$i %in% 2:length(mtch_sub)){
            
            # dercrease index
            v$i <- v$i - 1
            
            # chose new match
            x <- mtch_sub[v$i]
            
            # change to new text
            observe(output$t_number <- renderText(HTML(paste0("<b>", v$i, "/", length(mtch_sub), " ", x$target_name, "</b>"))))
            observe(output$t_formula <- renderText(HTML(paste0("Formula: ", x$target_formula, " InchiKey: ", x$target_inchikey))))
            observe(output$t_query <- renderText(HTML(paste0("Query: ", round(x$precursorMz,5), " @ ", round(x$rtime/60,2)))))
            observe(output$t_library <- renderText(HTML(paste0("Library: ", round(x$target_precursorMz,4), " @ ", round(x$target_rtime/60,2)))))
            observe(output$t_score <- renderText(HTML(paste0("Score: ", round(x$score, 3), " Reverse: ", round(x$reverse_score, 3)))))
            
            # render new plotly plot
            output$plot <-  renderPlotly(plotly_headtail(x))
            
            # change preselection in input table to new value => works, but somehow slows down the process
            # updateSelectInput(session = getDefaultReactiveDomain(), inputId = 'selection', selected=v$i)
            
            # change the button to the previous selected verification value
            updateRadioButtons(inputId="veri",
                 choices = list("True match" = T, "False match" = F), 
                 selected=r$veri_logical[v$i])
        }
    })
    
    # selection in list
    observeEvent(input$selection, {
            # change index to selection
            v$i <- as.numeric(input$selection)
            
            # chose new match
            x <- mtch_sub[v$i]
            
            # change to new text
            observe(output$t_number <- renderText(HTML(paste0("<b>", v$i, "/", length(mtch_sub), " ", x$target_name, "</b>"))))
            observe(output$t_formula <- renderText(HTML(paste0("Formula: ", x$target_formula, " InchiKey: ", x$target_inchikey))))
            observe(output$t_query <- renderText(HTML(paste0("Query: ", round(x$precursorMz,5), " @ ", round(x$rtime/60,2)))))
            observe(output$t_library <- renderText(HTML(paste0("Library: ", round(x$target_precursorMz,4), " @ ", round(x$target_rtime/60,2)))))
            observe(output$t_score <- renderText(HTML(paste0("Score: ", round(x$score, 3), " Reverse: ", round(x$reverse_score, 3)))))
            
            # render new plotly plot
            output$plot <-  renderPlotly(plotly_headtail(x))
            
            # change the button to the previous selected verification value
            updateRadioButtons(inputId="veri",
                               choices = list("True match" = T, "False match" = F), 
                               selected=r$veri_logical[v$i])
    })
    
    # Define what to do if click on radio button
    observeEvent(input$veri,{
        # store value in logical vector
        r$veri_logical[v$i] <- input$veri
        # update the names of the list
        names(l_choices) <- paste(seq(1, length(mtch_sub)),mtch_sub$target_name, "-", r$veri_logical)
        # update the list
        updateSelectInput(session = getDefaultReactiveDomain(), choices = l_choices, inputId = 'selection', selected=v$i)
    })
    
    #Store Result
    observeEvent(input$b_store,{
        # subset the match object
        mtch_subII <- mtch_sub[which(r$veri_logical == TRUE)]
        # store as RDS with changed name
        saveRDS(mtch_subII, paste0(str_replace(filename, ".rds$", "_verified.rds")))
    })
}
#######################################################################
# Run the application 
shinyApp(ui = ui, server = server)
