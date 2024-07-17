library(tidyverse)
library(shiny)
library(gridExtra)

final <- read_tsv("app/hgt_introns.tsv", col_names = T)


# Define UI
ui <- fluidPage(
  titlePanel("Density Plot and Box Plot of Exon-Exon Junctions - scaled positions"),
  
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
                  value = c(min(final$exons), max(final$exons)))
    ),
    
    mainPanel(
      plotOutput("combinedPlot")
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
    
    list(data1 = data1, data2 = data2)
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
      xlim(-1, 2) +  # Set the x-axis range
      scale_color_manual(values = c("blue", "red"))  # Match colors to geom_density fills
    
    # Box plots with jitter
    box_plot <- ggplot() +
      geom_boxplot(data = filtered_data()$data1, aes(x = scaled_pos, y = "Set 1"), fill = "blue", alpha = 0.1) +
      geom_jitter(data = filtered_data()$data1, aes(x = scaled_pos, y = "Set 1"), color = "blue", alpha = 0.1, height = 0.15) +
      geom_boxplot(data = filtered_data()$data2, aes(x = scaled_pos, y = "Set 2"), fill = "red", alpha = 0.1) +
      geom_jitter(data = filtered_data()$data2, aes(x = scaled_pos, y = "Set 2"), color = "red", alpha = 0.1, height = 0.15) +
      theme_minimal() +
      labs(x = "Scaled Gene Structure", y = "") +
      xlim(-1, 2) +  # Set the x-axis range
      theme(axis.title.y = element_blank(), axis.text.y = element_text(color = c("blue", "red")))
    
    # Combine the two plots vertically
    grid.arrange(density_plot, box_plot, ncol = 1, heights = c(2/3, 1/3))
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
