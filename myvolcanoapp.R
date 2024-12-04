#
# This is my pet project: a Shiny app that renders custom, dynamic volcano plots
# Date Created: December 4 2024
#

library(shiny)
library(ggplot2)
library(ggrepel)

make_volcano = function( 
    de, 
    log_fc_cutoff = 0.05, 
    p_val_cutoff = 0.01,  
    graph_title = 'Add a Title', 
    overlap_metric = 17, 
    label_genes = TRUE, 
    label_specific_genes = FALSE,
    genes_to_label = c(),
    upcolor = 'Yellow', 
    downcolor = 'Blue', 
    midcolor = 'Gray', 
    labelcolor = 'Black', 
    remove_genes = c(),
    remove_labels = c(),
    overlap_metric_secondary = Inf, 
    visible_cutoffs = TRUE) {
  
  if(!("gene" %in% colnames(de))) { warning("Differential Expression Dataframe Must Have `gene` Column!")}
  
  de = subset(de,!(gene %in% remove_genes))
  
  up_de <- subset(de, avg_log2FC > 0)
  down_de <- subset(de, avg_log2FC < 0)

  up_quant <- log_fc_cutoff
  down_quant <- -1 * log_fc_cutoff
  
  de$reg = ""
  de$reg[de$p_val_adj < p_val_cutoff & de$avg_log2FC > up_quant & de$avg_log2FC > 0] <- "UP"
  de$reg[de$p_val_adj < p_val_cutoff & de$avg_log2FC < down_quant & de$avg_log2FC < 0] <- "DOWN"
  de$name = de$gene
  de$name[de$reg == ""] <- ""
  de$name[de$name %in% remove_labels] <- ""
  if(label_specific_genes) {de <- de %>% mutate(name_specific = ifelse(gene  %in% genes_to_label & reg != "",gene, ""))}    
  
  
  plot = ggplot(data=de, aes(x=avg_log2FC, y=-log10(p_val_adj), col=reg, label=name)) + 
    geom_point(color = 'black', size = 2.5) + 
    theme_minimal() +
    scale_color_manual(breaks = c("DOWN", "", "UP"),values=c(downcolor, midcolor, upcolor)) + 
    geom_point(size = 1) + 
    ggtitle(graph_title) +
    theme(panel.grid = element_blank(), legend.position = "none") 
  
  if (label_genes) {plot = plot + ggrepel::geom_text_repel(aes(label = name), color = labelcolor, max.overlaps = overlap_metric, nudge_y = 1)}
  if (label_specific_genes) {plot = plot + geom_text_repel(aes(label = name_specific), color = labelcolor, max.overlaps = overlap_metric_secondary)}
  if (visible_cutoffs) {plot = plot + geom_vline(xintercept = up_quant, linetype=3) + geom_vline(xintercept = down_quant, linetype = 3) + geom_hline(yintercept = -log10(p_val_cutoff), linetype = 3)}
  
  return(plot)
}

# Define UI for application that draws a histogram
ui <- fluidPage(
  titlePanel("Volcano Plot"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("upload", "Upload Differential Expression Results (.csv file)", accept = c(".csv")),
      textInput("title", "Enter Title", value = "My Volcano Plot"),
      h4("Labels"),
      actionButton("labels","Show/Hide Labels"),
      numericInput("overlap","Overlap Metric", value = 10),
      h4("Thresholds"),
      numericInput("logfc", "Log Fold Change Threshold", value = 0.5),
      numericInput("pval", "FDR P Value Threshold", value = 0.01),
      h4("Colors"),
      textInput("upcolor","Upregulated Genes", value = "Yellow"),
      textInput("downcolor","DownregulatedGenes",value = "Blue"),
      textInput("midcolor","Non-significant Genes", value = "Grey")
      
    ),
    mainPanel(
      plotOutput("plot")
    )
  )
)

server <- function(input, output, session) {
  data <- reactive({
    req(input$upload)
    ext <- tools::file_ext(input$upload$name)
    switch(ext,
           csv = vroom::vroom(input$upload$datapath, delim = ","),
           validate("Invalid file; Please upload a .csv file")
    )
  })
  title <- reactive(input$title)
  
  logfc <- reactive(input$logfc)
  pval <- reactive(input$pval)
  
  overlap <- reactive(input$overlap)
  
  upcolor <- reactive(input$upcolor)
  downcolor <- reactive(input$downcolor)
  midcolor <- reactive(input$midcolor)
  
  
  # Toggle switch for showing labels
  show_labels <- reactiveVal(TRUE)
  observeEvent(input$labels, {
    show_labels(!show_labels())
  })
  
  output$plot <- renderPlot({
    make_volcano(data(), 
                 graph_title = title(),
                 overlap_metric = overlap(),
                 label_genes = show_labels(),
                 log_fc_cutoff = logfc(), 
                 p_val_cutoff = pval(),
                 upcolor = upcolor(), 
                 downcolor = downcolor(), 
                 midcolor = midcolor())
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
