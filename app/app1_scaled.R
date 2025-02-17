library(tidyverse)
library(shiny)
library(gridExtra)
library(DT)
library(cowplot)

final <- read_tsv("app/hgt_introns.tsv", col_names = TRUE)

# Define UI
ui <- fluidPage(
  titlePanel("Density Plot and Box Plot of Exon-Exon Junctions - Scaled gene structure"),
  
  sidebarLayout(
    sidebarPanel(
      # First set of input controls
      checkboxGroupInput("rank1", "Rank Set 1", 
                         choices = unique(final$rank),
                         selected = unique(final$rank)),
      sliderInput("cds_size1", "CDS Size Set 1", 
                  min = min(final$cds_size), max = 30000, 
                  value = c(min(final$cds_size), 30000)),
      sliderInput("cDNA_size1", "cDNA Size Set 1",
                  min = min(final$cDNA_size), max = max(final$cDNA_size),
                  value = c(min(final$cDNA_size), max(final$cDNA_size))),
      sliderInput("exons1", "Exons Set 1", 
                  min = min(final$exons), max = max(final$exons), 
                  value = c(min(final$exons), max(final$exons))),
      
      # Second set of input controls
      checkboxGroupInput("rank2", "Rank Set 2", 
                         choices = unique(final$rank),
                         selected = unique(final$rank)),
      sliderInput("cds_size2", "CDS Size Set 2", 
                  min = min(final$cds_size), max = 30000, 
                  value = c(min(final$cds_size), 30000)),
      sliderInput("cDNA_size2", "cDNA Size Set 2",
                  min = min(final$cDNA_size), max = max(final$cDNA_size),
                  value = c(min(final$cDNA_size), max(final$cDNA_size))),
      sliderInput("exons2", "Exons Set 2", 
                  min = min(final$exons), max = max(final$exons), 
                  value = c(min(final$exons), max(final$exons))),
      downloadButton("downloadData1", "Download Table for Set 1"),
      downloadButton("downloadData2", "Download Table for Set 2")
    ),
    
    mainPanel(
      plotOutput("combinedPlot"),
      DTOutput("table1"),
      DTOutput("table2")
    )
  )
)

# Define server logic
server <- function(input, output) {
  
  filtered_data <- reactive({
    # Filter data for first set of parameters
    data1 <- final %>%
      filter(rank %in% input$rank1) %>%
      filter(cds_size >= input$cds_size1[1] & cds_size <= input$cds_size1[2]) %>%
      filter(cDNA_size >= input$cDNA_size1[1] & cDNA_size <= input$cDNA_size1[2]) %>%
      filter(exons >= input$exons1[1] & exons <= input$exons1[2])
    
    # Filter data for second set of parameters
    data2 <- final %>%
      filter(rank %in% input$rank2) %>%
      filter(cds_size >= input$cds_size2[1] & cds_size <= input$cds_size2[2]) %>%
      filter(cDNA_size >= input$cDNA_size2[1] & cDNA_size <= input$cDNA_size2[2]) %>%
      filter(exons >= input$exons2[1] & exons <= input$exons2[2])
    
    # Create distinct datasets for boxplot
    data3 <- data1 %>%
      distinct(gene_id, .keep_all = TRUE)
    
    data4 <- data2 %>%
      distinct(gene_id, .keep_all = TRUE)
    
    list(data1 = data1, data2 = data2, data3 = data3, data4 = data4)
  })
  
  output$combinedPlot <- renderPlot({
    data1_count <- nrow(filtered_data()$data1)
    data2_count <- nrow(filtered_data()$data2)
    
    # Density plot
    density_plot <- ggplot() +
      geom_density(data = filtered_data()$data1, aes(x = scaled_pos, color = "Set 1"), fill = "blue", alpha = 0.1) +
      geom_density(data = filtered_data()$data2, aes(x = scaled_pos, color = "Set 2"), fill = "red", alpha = 0.1) +
      theme_minimal() +
      labs(title = "Density Plot of Exon-Exon Junctions", x = "Scaled Gene Structure", y = "Density") +
      xlim(-1, 2) + # Set the x-axis range
      scale_color_manual(values = c("blue", "red")) # Match colors to geom_density fills
    
    # Box plots with jitter
    box_plot <- ggplot() +
      geom_boxplot(data = filtered_data()$data1, aes(x = scaled_pos, y = "Set 1"), fill = "blue", alpha = 0.1) +
      geom_jitter(data = filtered_data()$data1, aes(x = scaled_pos, y = "Set 1"), color = "blue", alpha = 0.1, height = 0.15) +
      geom_boxplot(data = filtered_data()$data2, aes(x = scaled_pos, y = "Set 2"), fill = "red", alpha = 0.1) +
      geom_jitter(data = filtered_data()$data2, aes(x = scaled_pos, y = "Set 2"), color = "red", alpha = 0.1, height = 0.15) +
      theme_minimal() +
      labs(x = "Scaled Gene Structure", y = "") +
      xlim(-1, 2) + # Set the x-axis range
      theme(axis.title.y = element_blank(), axis.text.y = element_text(color = c("blue", "red")))
    
    # Boxplot for expression data
    expression_plot <- ggplot() +
      geom_boxplot(data = filtered_data()$data3, aes(x = "Set 1", y = log2_expression), fill = "blue", alpha = 0.5) +
      geom_boxplot(data = filtered_data()$data4, aes(x = "Set 2", y = log2_expression), fill = "red", alpha = 0.5) +
      theme_minimal() +
      labs(title = "Expression Data", x = "", y = "log2(Expression)")
    
    # Combine the plots
    combined_density_box <- plot_grid(density_plot, box_plot, ncol = 1, align = 'v', rel_heights = c(2/3, 1/3))
    plot_grid(expression_plot, combined_density_box, ncol = 2, rel_widths = c(1/3, 2/3))
  })
  
  output$table1 <- renderDT({
    datatable(filtered_data()$data1, options = list(pageLength = 5))
  })
  
  output$table2 <- renderDT({
    datatable(filtered_data()$data2, options = list(pageLength = 5))
  })
  
  output$downloadData1 <- downloadHandler(
    filename = function() {
      paste("filtered_data_set1-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(filtered_data()$data1, file, row.names = FALSE)
    }
  )
  
  output$downloadData2 <- downloadHandler(
    filename = function() {
      paste("filtered_data_set2-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(filtered_data()$data2, file, row.names = FALSE)
    }
  )
}

# Run the application 
shinyApp(ui = ui, server = server)

