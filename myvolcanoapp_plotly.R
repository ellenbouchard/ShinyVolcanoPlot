#
# This is my pet project: a Shiny app that renders custom, dynamic volcano plots
# Date Created: December 4 2024
#

library(shiny)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(plotly)

### FUNCTIONS
make_volcano_plotly = function( 
    de, 
    log_fc_cutoff = 0.05, 
    p_val_cutoff = 0.01,  
    graph_title = 'Add a Title', 
    overlap_metric = 17, 
    label_genes = TRUE, 
    genes_to_label = c(),
    upcolor = 'Yellow', 
    downcolor = 'Blue', 
    midcolor = 'Gray', 
    labelcolor = 'Black', 
    pointsize = 5,
    labelsize = 3,
    remove_genes = c(),
    remove_labels = c(),
    overlap_metric_secondary = Inf, 
    visible_cutoffs = TRUE,
    show_outlines = TRUE,
    height = 500,
    width = 400) {
  
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
 # de$name[de$reg == ""] <- ""
  #de$name[de$name %in% union(remove_labels, genes_to_label)] <- ""
  if(length(genes_to_label) > 0) {de <- de %>% mutate(name_specific = ifelse(gene  %in% genes_to_label & reg != "",gene, ""))}    
  
  #plot = ggplot(data=de, aes(x=avg_log2FC, y=neg_log10_pval, col=reg, label=name)) + 
  #  geom_point(color = 'black', size = 2.5) + 
  #  theme_minimal() +
  #  scale_color_manual(breaks = c("DOWN", "", "UP"),values=c(downcolor, midcolor, upcolor)) + 
  #  geom_point(size = 1) + 
  #  ggtitle(graph_title) +
  #  theme(panel.grid = element_blank(), legend.position = "none") 
  
  #if (label_genes) {plot = plot + ggrepel::geom_text_repel(aes(label = name), color = labelcolor, size = labelsize, max.overlaps = overlap_metric, nudge_y = 1)}
  #if (length(genes_to_label) > 0) {plot = plot + geom_text_repel(aes(label = name_specific), color = labelcolor, size = labelsize, max.overlaps = overlap_metric_secondary)}
  #if (visible_cutoffs) {plot = plot + geom_vline(xintercept = up_quant, linetype=3) + geom_vline(xintercept = down_quant, linetype = 3) + geom_hline(yintercept = -log10(p_val_cutoff), linetype = 3)}
  
  # Define colors
  colors_pal <- c(downcolor, midcolor, upcolor)
  colors_pal <- setNames(colors_pal, c("DOWN", "", "UP")) 
  
  # Define threshold lines
  # Horizontal line
  hline <- function(y = 0, color = "black") {
    list(
      type = "line",
      x0 = 0,
      x1 = 1,
      xref = "paper",
      y0 = y,
      y1 = y,
      line = list(color = color, dash = "dot", width = 1)
    )
  }
  
  # Vertical line
  vline <- function(x = 1, color = "black") {
    list(
      type = "line",
      y0 = 0,
      y1 = 1,
      yref = "paper",
      x0 = x,
      x1 = x,
      line = list(color = color, dash = "dot", width = 1)
    )
  }
  
  plot <- plot_ly(data = de, 
                  x = ~avg_log2FC, 
                  y = ~neg_log10_pval,
                  color = ~reg, 
                  colors = colors_pal,
                  text = ~gene)
  # Add black out lines if show_outlines is TRUE
  if(show_outlines) {
    plot <- plot %>% 
      add_trace(
        x = ~avg_log2FC, 
        y = ~neg_log10_pval, 
        marker = list(
          size = pointsize + 2, 
          color = 'black', 
          line = list(width = 0)
        ),
        type = 'scatter',
        mode = 'markers',
        showlegend = FALSE,
        hoverinfo = "none"
      )
  } # Add markers
  plot <- plot %>% 
     add_trace(
      x = ~avg_log2FC, 
      y = ~neg_log10_pval, 
      text = ~paste(name, " \n log2fc:", round(avg_log2FC, 2), " \n adjPval:", p_val_adj),
      marker = list(
        size = pointsize, 
        line = list(width = 0)),
      type = 'scatter',
      mode = 'markers',
      showlegend = FALSE,
      hovertemplate = "%{text}<extra></extra>"
    )
  
  # Add labels
  plot <- plot %>% add_annotations(
    x = de$avg_log2FC,
    y = de$neg_log10_pval,
    text = de$gene,
    xref = "x",
    yref = "y",
    showarrow = TRUE,
    arrowhead = 0,
    arrowsize = 0,
    arrowwidth = 1,
    ax = 20,
    ay = -20,
    font = list(color = 'black',
                size = 10),
    showlegend = FALSE
  )

  
  plot <- plot %>% layout(title = graph_title, 
                          height = height,
                          width = width,
                          xaxis = list(
                            showgrid = FALSE,
                            zeroline = FALSE 
                          ),
                          yaxis = list(
                            showgrid = FALSE,
                            zeroline = FALSE  
                          ))
  if(visible_cutoffs) {
    plot <- plot %>% layout(shapes = list(hline(-log10(p_val_cutoff)), vline(up_quant), vline(down_quant)))
  }
  
    
  
  return(plot)
}
# Function for determining if input color is valid
is_valid_color <- function(color) {
  tryCatch({
    grDevices::col2rgb(color)
    TRUE
  }, error = function(e) {
    FALSE
  })
}

## SHINY APP

# Define UI for application
ui <- fluidPage(
  # Enable ShinyFeedback
  shinyFeedback::useShinyFeedback(),
  
  titlePanel("Volcano Plot"),
  
  sidebarLayout(
    sidebarPanel(
      tabsetPanel(
        tabPanel("General", 
                 fileInput("upload", "Upload Differential Expression Results (.csv file)", accept = c(".csv")),
                 textInput("title", "Enter Title", value = "My Volcano Plot"),
                 h4("Thresholds"),
                 numericInput("logfc", "Log Fold Change Threshold", value = 0.5),
                 numericInput("pval", "FDR P Value Threshold", value = 0.01),
                 actionButton("thresholds","Show/Hide Thresholds"),
                 h4("Plot Size"),
                 sliderInput("height", "Height", min = 100, max = 1000, value = 500),
                 sliderInput("width", "Width", min = 100, max = 1000, value = 500)
        ),
        tabPanel("Points",
                 h4("Colors"),
                 textInput("upcolor","Upregulated Genes", value = "Yellow"),
                 textInput("downcolor","DownregulatedGenes",value = "Blue"),
                 textInput("midcolor","Non-significant Genes", value = "Grey"),
                 actionButton("outlines", "Show/Hide Outlines"),
                 sliderInput("pointsize", "Point Size", min = 1, max = 15, value = 5)
        ),
        tabPanel("Labels", 
                 actionButton("labels","Show/Hide Default Labels"),
                 numericInput("overlap","Default Label Density", value = 10),
                 h4("Add Custom Labels by Entering Gene Names"),
                 textAreaInput("gene_input", "Enter Gene Names (comma or newline separated)",
                               height = "100px"),
                 numericInput("labelsize", "Label Font Size", value = 4),
                 textInput("labelcolor", "Label Font Color", value = "Black")
        ),
        tabPanel("Download",
                 h4("Download Your Plot"),
                 selectInput("filetype", "Select File Type", choices = c("png","jpg","svg"), selected = "png"),
                 downloadButton("download")
        ),
      )
    ),
    mainPanel(
      # Add text instructions above plot
      fluidRow(
        div(
          style = "padding: 10px; background-color: #f8f9fa; border: 1px solid #dee2e6; border-radius: 5px; margin-bottom: 10px;",
          HTML("
            <ul>
              <li>Hover over a point to see gene info.</li>
            </ul>
          ")
        )
      ),
      fluidRow(
        plotlyOutput("plot")#,
        #tableOutput("data")
      )
    )
  )
)

server <- function(input, output, session) {
  # Process the input file
  data <- reactive({
    req(input$upload)
    ext <- tools::file_ext(input$upload$name)
    # Warning if the file is not a .csv
    # This is probably not necessary since the app only accepts .csv files anyway
    switch(ext,
           csv = vroom::vroom(input$upload$datapath, delim = ","),
           validate("Invalid file; Please upload a .csv file")
    )
    # Read the file into a dataframe
    df <- read.csv(input$upload$datapath)
    
    # Check for required columns
    required_cols <- c("gene", "avg_log2FC", "p_val_adj")
    missing_cols <- setdiff(required_cols, colnames(df))
    
    # Issue warning if required columns are missing
    shinyFeedback::feedbackWarning(
      inputId = "upload",
      show = length(missing_cols) > 0,
      text = paste("Missing required columns:", paste(missing_cols, collapse = ", "))
    )
    
    # Stop processing if required columns are missing
    validate(need(length(missing_cols) == 0, ""))
    
    # Create a separate column for -log10(p_val_adj)
    # This is necessary for enabling the hover and click functions on the plot
    df$neg_log10_pval <- -log10(df$p_val_adj)
    return(df)
  })
  
  # Define input variables as reactives
  title <- reactive(input$title)
  logfc <- reactive(input$logfc)
  pval <- reactive(input$pval)
  overlap <- reactive(input$overlap)
  labelsize <- reactive(input$labelsize)
  pointsize <- reactive(input$pointsize)
  height <- reactive(input$height)
  width <- reactive(input$width)
  
  # For colors, issue warning if color is not valid
  upcolor <- reactive({
    validcolor <- is_valid_color(input$upcolor)
    shinyFeedback::feedbackWarning("upcolor", !validcolor, "Please select a valid color or HEX code")
    req(validcolor)
    input$upcolor
    })
  downcolor <- reactive({
    validcolor <- is_valid_color(input$downcolor)
    shinyFeedback::feedbackWarning("downcolor", !validcolor, "Please select a valid color or HEX code")
    req(validcolor)
    input$downcolor
  })
  midcolor <- reactive({
    validcolor <- is_valid_color(input$midcolor)
    shinyFeedback::feedbackWarning("midcolor", !validcolor, "Please select a valid color or HEX code")
    req(validcolor)
    input$midcolor
    })
  labelcolor <- reactive({
    validcolor <- is_valid_color(input$labelcolor)
    shinyFeedback::feedbackWarning("labelcolor", !validcolor, "Please select a valid color or HEX code")
    req(validcolor)
    input$labelcolor
  })
  
  
  # Toggle switch for showing default labels
  show_labels <- reactiveVal(TRUE)
  observeEvent(input$labels, {
    show_labels(!show_labels())
  })
  
  # Toggle switch for showing thresholds
  show_cutoffs <- reactiveVal(TRUE)
  observeEvent(input$thresholds, {
    show_cutoffs(!show_cutoffs())
  })
  
  # Toggle switch for showing outlines
  show_outlines <- reactiveVal(TRUE)
  observeEvent(input$outlines, {
    show_outlines(!show_outlines())
  })
  
  # Make reactive lists to store specific genes to label and specific genes to remove
  clicked_on_genes <- reactiveVal(c())
  genes_to_label <- reactiveVal(c())
  remove_labels <- reactiveVal(c())
  
  # Process user input to handle comma, space, and newline separation
  observeEvent(input$gene_input, {
    # Split input by comma, space, and newline, and remove empty strings
    input_text <- gsub("\n", " ", input$gene_input)
    input_genes <- unlist(strsplit(input$gene_input, "[,\n]+"))
    input_genes <- trimws(input_genes)
    input_genes <- input_genes[input_genes != ""]
    
    # Update genes_to_label dynamically
    # Genes to label should include all genes in the text input, 
    # As well as all genes that are currently clicked on
    all_genes <- unique(c(input_genes, clicked_on_genes()))
    genes_to_label(all_genes)
  })
  
  # Observe clicks on the plot
  observeEvent(input$plot_click, {
    click_info <- nearPoints(data(), input$plot_click, xvar = "avg_log2FC", yvar = "neg_log10_pval")
    if (nrow(click_info) > 0) {
      clicked_gene <- click_info$gene[1]
      current_clicked_genes <- clicked_on_genes() # Get list of genes that are already clicked on
      
      # Get the current list of genes in genes_to_label (both from input and clicks)
      current_genes_to_label <- genes_to_label()
      input_text <- gsub("\n", " ", input$gene_input)
      input_genes <- unlist(strsplit(input$gene_input, "[,\n]+"))
      input_genes <- trimws(input_genes)
      input_genes <- input_genes[input_genes != ""]
      
      # If the gene is not already in the list of clicked "on" genes OR in the input genes, add it
      # Note: you have to make sure the gene isn't already in the text input to avoid wonkiness
      if (!(clicked_gene %in% current_clicked_genes || clicked_gene %in% input_genes)) {
        clicked_on_genes(unique(c(current_clicked_genes, clicked_gene)))
        # If gene is already in current_clicked_genes, but not in input_genes, remove it
        # Also update remove_genes
      } else if (clicked_gene %in% current_clicked_genes && !(clicked_gene %in% input_genes)) { 
        clicked_on_genes(setdiff(current_clicked_genes, clicked_gene))
        remove_labels(c(remove_labels, clicked_gene))
      }
      
      # Update genes_to_label with input genes and clicked "on" genes
      all_genes <- unique(c(input_genes, clicked_on_genes()))
      genes_to_label(all_genes)

      
    }
  })
  
  
  # Render the plot
  output$plot <- renderPlotly(
    {
    p <- make_volcano_plotly(data(), 
                 graph_title = title(),
                 overlap_metric = overlap(),
                 label_genes = show_labels(),
                 log_fc_cutoff = logfc(), 
                 p_val_cutoff = pval(),
                 upcolor = upcolor(), 
                 downcolor = downcolor(), 
                 midcolor = midcolor(),
                 genes_to_label = genes_to_label(),
                 remove_labels = remove_labels(),
                 labelcolor = labelcolor(),
                 labelsize = labelsize(),
                 pointsize = pointsize(),
                 visible_cutoffs = show_cutoffs(),
                 show_outlines = show_outlines(),
                 height = height(),
                 width = width())
    ggplotly(p) %>% config(editable = TRUE)
  })
  
  # Render the table
  output$data <- renderTable({
    selected_data <- nearPoints(data(), input$plot_hover, xvar = "avg_log2FC", yvar = "neg_log10_pval")
    # Format p values in scientific notation
    selected_data$p_val_adj <- format(selected_data$p_val_adj, scientific = TRUE, digits = 3)
    selected_data$p_val <- format(selected_data$p_val, scientific = TRUE, digits = 3)
    selected_data
  })
  
  # Render the figure for download
  output$download <- downloadHandler(
    filename = function() {
      paste("plot", Sys.Date(), ".", input$filetype, sep = "")
    },
    content = function(file) {
      ggsave(file, plot = make_volcano(data(), 
                                       graph_title = title(),
                                       overlap_metric = overlap(),
                                       label_genes = show_labels(),
                                       log_fc_cutoff = logfc(), 
                                       p_val_cutoff = pval(),
                                       upcolor = upcolor(), 
                                       downcolor = downcolor(), 
                                       midcolor = midcolor(),
                                       genes_to_label = genes_to_label(),
                                       remove_labels = remove_labels(),
                                       labelcolor = labelcolor(),
                                       labelsize = labelsize())
            , device = input$filetype, width = input$width / 100, height = input$height / 100)
    }
  )
}

# Run the application 
shinyApp(ui = ui, server = server)
