
library(shiny)
library(shinyWidgets)
library(DT)
library(tidyverse)
library(plotly)


TemplateExample <- read_csv("TemplateExample.csv")


ui <- bootstrapPage(

    #load in the font from google 
    tags$link(href='https://fonts.googleapis.com/css?family=Montserrat', rel='stylesheet'),
    
    #Lots of CSS
    tags$head(tags$style(HTML('
    
                                body {
                               font-family:"Montserrat";
                                }
                              .navbar-bg {
                               background-color: #003776 !important;
                              }
                              .nav > li > a:hover, .nav > li > a:focus {
                              background-color: #FFFFFF;
                                }
                              .bg-grey {
                                background-color: #d6d6d6;
                              }
                              .bg-dark-grey {
                                background-color: #696868;
                                color: white;
                              }
                              hr{
                                height: 3px;
                                background-color: #014093;
                                border: none;
                                width:10%;
                              }
                              
                              .sidenav {
                                height: 100%; /* Full-height: remove this if you want "auto" height */
                                width: 230px; /* Set the width of the sidebar */
                                position: fixed; /* Fixed Sidebar (stay in place on scroll) */
                                z-index: 1; /* Stay on top */
                                top: 0; /* Stay at the top */
                                left: 0;
                                background-color: #111; /* Black */
                                overflow-x: hidden; /* Disable horizontal scroll */
                                padding-top: 20px;
                              }

                              /* The navigation menu links */
                              .sidenav a {
                                padding: 6px 8px 6px 16px;
                                text-decoration: none;
                                font-size: 20px;
                                color: #818181;
                                display: block;
                              }

                              /* When you mouse over the navigation links, change their color */
                              .sidenav a:hover {
                                color: #f1f1f1;
                              }


                              /* On smaller screens, where height is less than 450px, change the style of the sidebar (less padding and a smaller font size) */
                              @media screen and (max-height: 450px) {
                                .sidenav {padding-top: 15px;}
                                .sidenav a {font-size: 18px;}
                              }
                              
                              
                              '))),
    tags$style(type="text/css",
               ".shiny-output-error { visibility: hidden; }",
               ".shiny-output-error:before { visibility: hidden; }"
    ),
    
    div(class= "fluid-page",
        div(class = "row text-center",
            div(class = "col-xl-12",style = "margin-left: 180px; margin-right: 180px;",
                br(),
                strong(h1('Multiplex Analysis')),
                hr(),
                uiOutput("controls_ui"),
                hr(),
                br(),
                br(),
                radioGroupButtons(
                    inputId = "display_type",
                    label = "",
                    #This is all the types of pages
                    choices = c("Base Graph", "Target", "Density"),
                    status = "info",
                    justified = TRUE,
                    checkIcon = list(
                        yes = icon("ok", 
                                   lib = "glyphicon"))
                ),
                br(),
                br(),
                #This is the graph output
                column(12, align = "center",
                plotlyOutput("graph", width = "90%", height = 600),
                br(),
                br(),
                br(),
                br(),
                downloadButton('downloadPlot', 'Download Plot'),
                br(),
                br(),
                br(),
                h1("Difference Summaries"),
                hr(),
                h4("Base Level is cycles 5-10, End Level is cycles 50-end "),
                br(),
                br(),
                DTOutput("sum_dat"),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br(),
                br()
                )
                
                
            )
            
            
        )
    )
    
    
    
    )

server <- function(input, output) {

    #When the app is open enter information
    
    showModal(modalDialog(
        title = "Enter Multiplex Data",
        easyClose = FALSE,
        footer = NULL,
        size = "m",
        fileInput("csvFile", "Input Flouresence File", width = 1200),
        fileInput("csvFile2", "Input Plate File", width = 1200),
        downloadLink('downloadData', 'Download Template/Example Plate File'),
        br(),
        br(),
        actionButton("submit_files", label = "Done")
    ))
    
    
    #Get the raw data from upload
    rawData <- reactive({
        inFile <- input$csvFile
        if (is.null(inFile)) return(NULL)
        data <- read.csv(inFile$datapath, header = TRUE)
        data
    })
    
    #Get the template from upload
    rawData2 <- reactive({
        inFile <- input$csvFile2
        if (is.null(inFile)) return(NULL)
        data <- read.csv(inFile$datapath, header = TRUE)
        data
    })
    
    #Merge them bois
    full_data <- reactive({
        TitrationLevels <- rawData()
        TemplateExample <- rawData2()
        
        colnames(TitrationLevels)[1] <- "Cycle"
        
        titration_long <- TitrationLevels %>% 
            pivot_longer(-Cycle)
        
        full_data <- left_join(titration_long, TemplateExample, by = c("name" = "Sample"))
        
        colnames(full_data)[4] <- "Probe concentration"
        colnames(full_data)[6] <- "Target Concentration"
        
        
        full_data
    })
    
    output$sum_dat <- renderDT({
      full_data <- full_data()
      
      
      Name <- unique(full_data$name)
      base_level <- rep(0,length(Name))
      end_level <- rep(0,length(Name))
      target <- rep(0,length(Name))
      index <- 1
      for (i in Name) {
        name_dat <- full_data[full_data$name == i,]
        target[index] <- name_dat$Target[1]
        low <- name_dat$value[name_dat$Cycle >=5 & name_dat$Cycle <= 10]
        high <- name_dat$value[name_dat$Cycle >=50]
        base_level[index] <- round(mean(low),2)
        end_level[index] <- round(mean(high),2)
        index = index + 1
      }
      difference <- end_level - base_level
      sum_dat <- data.frame(Name, target, base_level, end_level, difference)
      datatable(sum_dat, filter = "top")
      
      
    })
    
    
                                
    #Download a template
    output$downloadData <- downloadHandler(
           filename = function() {
             paste('TemplateExample-', Sys.Date(), '.csv', sep='')
           },
           content = function(con) {
             write.csv(TemplateExample, con)
           }
         )
    
   #Remove modal when done 
    observeEvent(input$submit_files, {
        
        removeModal()
        
    })
    
    #Check that tables are working correctly 
    output$table <- renderTable({
        
        full_data_summary()
        
    })
    
    output$controls_ui <- renderUI({
        
        full_data <- full_data()
        
        tagList(
            column(12, align="center",
            pickerInput(
                inputId = "filter_target",
                label = "Filter Targets", 
                choices = as.character(unique(full_data$Target)),
                selected = as.character(unique(full_data$Target)),
                options = list(
                    `actions-box` = TRUE), 
                multiple = TRUE
            ),
            
            pickerInput(
                inputId = "filter_samples",
                label = "Filter Samples", 
                choices = unique(full_data$name),
                selected = unique(full_data$name),
                options = list(
                    `actions-box` = TRUE,
                    `live-search` = TRUE), 
                multiple = TRUE
            ),
            
            pickerInput(
                inputId = "filter_tc",
                label = "Filter Target Concentration", 
                choices = as.character(unique(full_data$`Target Concentration`)),
                selected = as.character(unique(full_data$`Target Concentration`)),
                options = list(
                    `actions-box` = TRUE,
                    `live-search` = TRUE), 
                multiple = TRUE
            ),
            
            pickerInput(
                inputId = "filter_pc",
                label = "Filter Probe Concentration", 
                choices = as.character(unique(full_data$`Probe concentration`)),
                selected = as.character(unique(full_data$`Probe concentration`)),
                options = list(
                    `actions-box` = TRUE,
                    `live-search` = TRUE), 
                multiple = TRUE
            ),
            
            sliderInput("filter_cycles", "Filter Cycles",
                        min = min(full_data$Cycle), 
                        max = max(full_data$Cycle), 
                        value =  c(min(full_data$Cycle),max(full_data$Cycle) ),
                        step = 1
            )
            
            
            
            )
            
        )
        
    })
    
    output$graph <- renderPlotly({
        
        
        
        
        full_data <- full_data()
        
        #filter Samples
        full_data <- full_data[full_data$name %in% input$filter_samples,]
        #Filter Targets
        full_data <- full_data[full_data$Target %in% input$filter_target,]
        #Filter Target Concentration
        full_data <- full_data[full_data$`Target Concentration` %in% input$filter_tc,]
        #Filter probe Concnetration
        full_data <- full_data[full_data$`Probe concentration` %in% input$filter_pc,]
        #Filter Cycles
        full_data <- full_data[full_data$Cycle <= input$filter_cycles[2],]
        full_data <- full_data[full_data$Cycle >= input$filter_cycles[1],]
        #Set names in correct order
        full_data$name <- factor(full_data$name, levels = unique(full_data$name))
        
        selection <- input$display_type
        
        
        if (selection == "Base Graph") {
            ggplotly(
            ggplot(data = full_data) +
                geom_line(aes(x = Cycle, y = value, color =  name, group = Target)) +
                theme_minimal() +
                ylab("")
            ) %>% 
                config(displayModeBar = F) 
        }
        else if (selection == "Target") {
            ggplotly(
                ggplot(data = full_data) +
                    geom_point(aes(x = Cycle, y = value, color =  Target)) +
                    theme_minimal() +
                    ylab("")
            ) %>% 
                config(displayModeBar = F) 
        }
        else {
            
            full_data_ends <- full_data[full_data$Cycle == max(full_data$Cycle),]
            
            full_data_summary <- full_data_ends %>% 
                group_by(Target) %>% 
                summarise(sd = sd(value), avg = mean(value))
            
            ggplotly(
            ggplot() +
                geom_density(data = full_data_ends, aes(x = value, fill = Target), alpha = 0.5) +
                scale_y_reverse() +
                coord_flip() +
                theme_minimal() +
                geom_point(data = full_data_ends, aes(x = value, y = -.001, color = Target)) +
                geom_text(data = full_data_summary, aes(x = avg, y = -.003, label = round(sd, 2))) +
                geom_text(aes(x = 900, y = -.003, label = "std dev")) +
                ylab("")
            )%>% 
                config(displayModeBar = F) 
            
        }
    })
    
    
    plotInput <- reactive({
      
      full_data <- full_data()
      
      #filter Samples
      full_data <- full_data[full_data$name %in% input$filter_samples,]
      #Filter Targets
      full_data <- full_data[full_data$Target %in% input$filter_target,]
      #Filter Target Concentration
      full_data <- full_data[full_data$`Target Concentration` %in% input$filter_tc,]
      #Filter probe Concnetration
      full_data <- full_data[full_data$`Probe concentration` %in% input$filter_pc,]
      #Filter Cycles
      full_data <- full_data[full_data$Cycle <= input$filter_cycles[2],]
      full_data <- full_data[full_data$Cycle >= input$filter_cycles[1],]
      #Order names
      full_data$name <- factor(full_data$name, levels = unique(full_data$name))
      
      selection <- input$display_type
      
      
      if (selection == "Base Graph") {
        
          ggplot(data = full_data) +
            geom_line(aes(x = Cycle, y = value, color =  name)) +
            theme_minimal() +
            ylab("")
        
      }
      else if (selection == "Target") {
        
          ggplot(data = full_data) +
            geom_point(aes(x = Cycle, y = value, color =  Target)) +
            theme_minimal() +
            ylab("")
        
      } 
      else {
        
        full_data_ends <- full_data[full_data$Cycle == max(full_data$Cycle),]
        
        full_data_summary <- full_data_ends %>% 
          group_by(Target) %>% 
          summarise(sd = sd(value), avg = mean(value))
        
        
          ggplot() +
            geom_density(data = full_data_ends, aes(x = value, fill = Target), alpha = 0.5) +
            scale_y_reverse() +
            coord_flip() +
            theme_minimal() +
            geom_point(data = full_data_ends, aes(x = value, y = -.001, color = Target)) +
            geom_text(data = full_data_summary, aes(x = avg, y = -.003, label = round(sd, 2))) +
            geom_text(aes(x = 900, y = -.003, label = "std dev")) +
            ylab("")
      
      }
      
    })
    
    
    output$downloadPlot <- downloadHandler(
      filename = function() { paste('plot', '.png', sep='') },
      content = function(file) {
        ggsave(file,plotInput())
      }
    )
    
}


shinyApp(ui = ui, server = server)
