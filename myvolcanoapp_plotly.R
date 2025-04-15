# 
# My Volcano: A Shiny App for custom Volcano Plot generation 
#
# Author: Ellen Bouchard
# 
# Date Created: December 4 2024
# 
# This is my pet project: a Shiny app that renders custom, dynamic, and interactive volcano plots with Plotly. 
# 
#

############# SETUP ################
library(shiny)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(plotly)
library(shinyjs)

############## FUNCTIONS #################

## Function: process_for_volcano
# This function is used to process data for plotting. 
# INPUT:
#   de = an R dataframe containing DE results
#   log_fc_cutoff = a numeric indicating log2 fold change threshold for display
#   p_val_cutoff = a numeric indicating FDR adjusted p value threshold for display
#   genes_to_label_1 = a vector of strings indicating which genes to add text labels to in the "first" label layer
#   genes_to_label_2 = a vector of strings indicating which genes to add text labels to in the "second" label layer
#   genes_to_label_click = a vector of strings indicating which genes to add text labels to that were chosen via direct click 

process_for_volcano = function(de,
                             log_fc_cutoff = 0.05,
                             p_val_cutoff = 0.01,
                             genes_to_label_1 = c(),
                             genes_to_label_2 = c(),
                             genes_to_label_click = c()
                             ) {
  
  # Set positive and negative log fold change cutoffs 
  up_cutoff <- log_fc_cutoff
  down_cutoff <- -1 * log_fc_cutoff
  
  # For adjusted p values of zero, change neg_log10_pval so that it is not infinite (change to machine minimum representable value)
  de$neg_log10_pval[de$p_val_adj == 0] <- -log10(.Machine$double.xmin)
  
  # Label genes as significantly "UP" or "DOWN" based on p val and logFC thresholds
  # do so by adding "reg" column that contains either "UP", "DOWN" or blanks for non-significant genes 
  de$reg = ""
  de$reg[de$p_val_adj < p_val_cutoff & de$avg_log2FC > up_cutoff & de$avg_log2FC > 0] <- "UP"
  de$reg[de$p_val_adj < p_val_cutoff & de$avg_log2FC < down_cutoff & de$avg_log2FC < 0] <- "DOWN"
  
  # Add first layer of labels by adding gene name to "name1" column 
  de$name1 = ""
  if(length(genes_to_label_1) > 0) {de <- de %>% mutate(name1 = ifelse(gene  %in% genes_to_label_1 & reg != "",gene, ""))}  
  
  # Add second layer of labels by adding gene name to "name2" column 
  de$name2 = ""
  if(length(genes_to_label_2) > 0) {de <- de %>% mutate(name2 = ifelse(gene  %in% genes_to_label_2 & reg != "",gene, ""))}  
  
  # Add third ("clicked") layer of labels, only if genes are not in first two layers of labels, by adding gene name to "nameclick" column 
  de$nameclick = ""
  if(length(genes_to_label_click) > 0) {
    de <- de %>% mutate(nameclick = ifelse(gene %in% genes_to_label_click & 
                                             reg != "" &
                                             !(gene %in% genes_to_label_1 | gene %in% genes_to_label_2), 
                                           gene, ""))
    }
  return(de)
}

# Function make_volcano_ploty: 
# Takes the output of process_for_volcano and returns volcano plot made with Plotly
# INPUT:
#   de = an R dataframe of DE results, previously processed by the process_for_volcano function
#   annot_de_1 = a dataframe that contains the avg_log2FC and neg_log10_pval values, as well as names, for 
#     genes annotated in layer 1
#   annot_de_2 = a dataframe that contains the avg_log2FC and neg_log10_pval values, as well as names, for 
#     genes annotated in layer 2
#   annot_de_click = a dataframe that contains the avg_log2FC and neg_log10_pval values, as well as names, for 
#     genes annotated via direct click 
#   log_fc_cutoff = a numeric indicating log2 fold change threshold for display
#   p_val_cutoff = a numeric indicating FDR adjusted p value threshold for display
#   graph_title = a string indicating graph title
#   upcolor = color to use for significantly upregulated genes
#   downcolor = color to use for significantly downregulated genes
#   midcolor = color to use for genes that are not significant
#   labelcolor1 = color for text of labels in annotation layer 1
#   labelcolor2 = color for text of labels in annotation layer 2
#   labelcolorclick = color for text of labels in annotation layer "cick" (labels chosen by direct click)
#   pointsize = numeric indicating size of points
#   labelsize1 = numeric indicating size of labels for annotation layer 1 
#   labelsize2 = numeric indicating size of labels for annotation layer 2
#   labelsizeclick = numeric indicating size of labels for annotation layer "click"
#   visiblecutoffs = boolean indicating whether to display dotted lines indicating significance thresholds
#   threshold_color = color for lines of significance thresholds 
#   show_outlines = boolean indicating whether points have distinct outlines
#   outline_color = color for outlines of points
#   height = numeric indicating height of graph,
#   width = numeric indicating width of graph 
#   space = numeric indicating amount of space between a labeled point and the line connecting the point to its label
#   download.filetype = indicates what file type the plot should download as
# OUTPUT:
#     plot = a Plotly object
make_volcano_plotly = function( 
    de, 
    annot_de_1,
    annot_de_2,
    annot_de_click,
    log_fc_cutoff = 0.05, 
    p_val_cutoff = 0.01,  
    graph_title = 'Add a Title', 
    upcolor = 'Yellow', 
    downcolor = 'Blue', 
    midcolor = 'Gray', 
    labelcolor1 = 'Black', 
    labelcolor2 = 'Red',
    labelcolorclick = "Black",
    pointsize = 5,
    labelsize1 = 10,
    labelsize2 = 10,
    labelsizeclick = 10,
    visible_cutoffs = TRUE,
    threshold_color = 'black',
    show_outlines = TRUE,
    outline_color = "black",
    height = 500,
    width = 400,
    space = 5,
    download.filetype = "png") {
  
  # Define positive and negative log2 fold change thresholds 
  up_cutoff <- log_fc_cutoff
  down_cutoff <- -1 * log_fc_cutoff

  # Define colors for downregulated, upregulated, and insignificant genes 
  colors_pal <- c(downcolor, midcolor, upcolor)
  colors_pal <- setNames(colors_pal, c("DOWN", "", "UP")) 
  
  # Define threshold lines
  # Horizontal line
  hline <- function(y = 0, color = threshold_color) {
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
  vline <- function(x = 1, color = threshold_color) {
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
  
  # Base code for plot
  plot <- plot_ly(data = de, 
                  height = height,
                  width = width,
                  x = ~avg_log2FC, 
                  y = ~neg_log10_pval,
                  color = ~reg, 
                  colors = colors_pal) %>%
            config(editable = TRUE, 
                   displaylogo = FALSE, 
                   modeBarButtonsToRemove = c(
                     'sendDataToCloud', 'autoScale2d', 'resetScale2d', 'toggleSpikelines',
                     'hoverClosestCartesian', 'hoverCompareCartesian',
                     'zoom2d','pan2d','select2d','lasso2d','zoomIn2d','zoomOut2d'
                   ),
                   toImageButtonOptions = list(format = download.filetype))
  
  # Add thresholds if toggled on
  if(visible_cutoffs) {
    plot <- plot %>% layout(shapes = list(hline(-log10(p_val_cutoff)), vline(up_cutoff), vline(down_cutoff)))
  }
  
  # Add annotations 1
  if(nrow(annot_de_1) > 0) {
    plot <- plot %>% add_annotations(
      x = annot_de_1$avg_log2FC,
      y = annot_de_1$neg_log10_pval,
      text = annot_de_1$name1,
      xref = "x",
      yref = "y",
      showarrow = TRUE,
      arrowhead = 0,
      arrowsize = 0,
      arrowwidth = 1,
      standoff = pointsize + (space - 2),
      ax = annot_de_1$ax,
      ay = annot_de_1$ay,
      font = list(color = labelcolor1,
                  size = labelsize1),
      showlegend = FALSE
    )
  }
  
  # Add annotations 2
  if(nrow(annot_de_2) > 0) {
    plot <- plot %>% add_annotations(
      x = annot_de_2$avg_log2FC,
      y = annot_de_2$neg_log10_pval,
      text = annot_de_2$name2,
      xref = "x",
      yref = "y",
      showarrow = TRUE,
      arrowhead = 0,
      arrowsize = 0,
      arrowwidth = 1,
      standoff = pointsize + (space - 2),
      ax = annot_de_2$ax,
      ay = annot_de_2$ay,
      font = list(color = labelcolor2,
                  size = labelsize2),
      showlegend = FALSE
    )
  }
  
  # Add clicked annotations 
  if(nrow(annot_de_click) > 0) {
    plot <- plot %>% add_annotations(
      x = annot_de_click$avg_log2FC,
      y = annot_de_click$neg_log10_pval,
      text = annot_de_click$nameclick,
      xref = "x",
      yref = "y",
      showarrow = TRUE,
      arrowhead = 0,
      arrowsize = 0,
      arrowwidth = 1,
      standoff = pointsize + (space - 2),
      ax = annot_de_click$ax,
      ay = annot_de_click$ay,
      font = list(color = labelcolorclick,
                  size = labelsizeclick),
      showlegend = FALSE
    )
  }

  # Add outlines to points if show_outlines is TRUE
  if(show_outlines) {
    plot <- plot %>% 
      add_trace(
        x = ~avg_log2FC, 
        y = ~neg_log10_pval, 
        marker = list(
          size = pointsize + 2, 
          color = outline_color, 
          line = list(width = 0)
        ),
        type = 'scatter',
        mode = 'markers',
        showlegend = FALSE,
        hoverinfo = "none"
      )
    
  } # Add markers / points
  plot <- plot %>% 
     add_trace(
      x = ~avg_log2FC, 
      y = ~neg_log10_pval, 
      text = ~paste(gene, " \n log2fc:", round(avg_log2FC, 2), " \n adjPval:", p_val_adj),
      marker = list(
        size = pointsize, 
        line = list(width = 0)),
      type = 'scatter',
      mode = 'markers',
      showlegend = FALSE,
      hovertemplate = "%{text}<extra></extra>"
    )

  # Finalize layout
  plot <- plot %>% layout(title = graph_title, 
                          xaxis = list(
                            title = 'Log2 Fold Change',
                            showgrid = FALSE,
                            zeroline = FALSE 
                          ),
                          yaxis = list(
                            title = '-log10 Adjusted P Value',
                            showgrid = FALSE,
                            zeroline = FALSE  
                          ),
                          dragmode = FALSE,
                          margin = list(
                            l = 70,
                            r = 70,
                            t = 70,
                            b = 70
                          ))
  return(plot)
}


# Function is_valid_color
# Determines if input color string is valid to use as a color for plotting 
# INPUT:
#   color = a string indicating what the user has input for a color
# OUTPUT:
#   boolean indicating whether the input string can be used as a color 
is_valid_color <- function(color) {
  tryCatch({
    grDevices::col2rgb(color)
    TRUE
  }, error = function(e) {
    FALSE
  })
}



############ SHINY APP #############

############ Define UI for application ###############
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
                 textInput("thresholdcolor", "Threshold Color", value = "Black"),
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
                 textInput("outlinecolor","Outlines", value = "Black"),
                 actionButton("outlines", "Show/Hide Outlines"),
                 sliderInput("pointsize", "Point Size", min = 1, max = 15, value = 5)
        ),
        tabPanel("Labels", 
                 h4("Edit Labels Added by Clicking"),
                 numericInput("labelsizeclick", "Label Font Size", value = 10),
                 textInput("labelcolorclick", "Label Font Color", value = "Black"),
                 actionButton("reset_clicked_labels", "Reset Labels"),
                 h4("Add Labels by Entering Gene Names"),
                 textAreaInput("gene_input_1", "Enter Gene Names (comma or newline separated)",
                               height = "100px"),
                 numericInput("labelsize1", "Label Font Size", value = 10),
                 textInput("labelcolor1", "Label Font Color", value = "Red"),
                 h4("Add Secondary Labels"),
                 textAreaInput("gene_input_2", "Enter Gene Names (comma or newline separated)",
                               height = "100px"),
                 numericInput("labelsize2", "Label Font Size", value = 10),
                 textInput("labelcolor2", "Label Font Color", value = "Blue"),
                 sliderInput("space", "Space between line and point", min = 0, max = 10, value = 5)
        ),
        tabPanel("Download",
                 h4("Click the camera icon on your plot to download"),
                 selectInput("filetype", "Select File Type", choices = c("png","svg"), selected = "png")
        )
      )
    ),
    mainPanel(
      # Add text instructions above plot
      fluidRow(
        div(
          style = "padding: 10px; background-color: #f8f9fa; border: 1px solid #dee2e6; border-radius: 5px; margin-bottom: 10px;",
          HTML("
            <ul>
              <li>Upload a .csv file with your differential expression results.</li>
              <li>Hover over a point to see gene info.</li>
              <li>Click on a point to add a label. Click again to remove.</li>
              <li>Enter gene names in the Labels tab to add labels. These can't be removed by clicking.</li>
              <li>Drag and drop labels to reposition.</li>
            </ul>
          ")
        )
      ),
      fluidRow(
        plotlyOutput("plot")
      )
    )
  )
)


################ SERVER FUNCTION ##################
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
    
    # Calculate -log10 pvalue and add as column. This is necessary for hover and click behaviors.
    df$neg_log10_pval <- -log10(df$p_val_adj)
    return(df)
  })
  
  # Define input variables as reactive values
  title <- debounce(reactive(input$title), millis = 500)
  logfc <- debounce(reactive(input$logfc), millis = 300)
  pval <- debounce(reactive(input$pval), millis = 300)
  labelsize1 <- debounce(reactive(input$labelsize1), millis = 300)
  labelsize2 <- debounce(reactive(input$labelsize2), millis = 300)
  labelsizeclick <- debounce(reactive(input$labelsizeclick), millis = 3000)
  pointsize <- debounce(reactive(input$pointsize), millis = 300)
  height <- debounce(reactive(input$height), millis = 300)
  width <- debounce(reactive(input$width), millis = 300)
  space <- debounce(reactive(input$space), millis = 300)
  filetype <- debounce(reactive(input$filetype), millis = 300)
  
  # For colors, issue warning if color is not valid
  upcolor <- reactive({
    validcolor <- is_valid_color(trimws(input$upcolor))
    shinyFeedback::feedbackWarning("upcolor", !validcolor, "Please select a valid color or HEX code")
    req(validcolor)
    input$upcolor
    })
  downcolor <- reactive({
    validcolor <- is_valid_color(trimws(input$downcolor))
    shinyFeedback::feedbackWarning("downcolor", !validcolor, "Please select a valid color or HEX code")
    req(validcolor)
    input$downcolor
  })
  midcolor <- reactive({
    validcolor <- is_valid_color(trimws(input$midcolor))
    shinyFeedback::feedbackWarning("midcolor", !validcolor, "Please select a valid color or HEX code")
    req(validcolor)
    input$midcolor
    })
  labelcolor1 <- reactive({
    validcolor <- is_valid_color(trimws(input$labelcolor1))
    shinyFeedback::feedbackWarning("labelcolor1", !validcolor, "Please select a valid color or HEX code")
    req(validcolor)
    input$labelcolor1
  })
  labelcolor2 <- reactive({
    validcolor <- is_valid_color(trimws(input$labelcolor2))
    shinyFeedback::feedbackWarning("labelcolor2", !validcolor, "Please select a valid color or HEX code")
    req(validcolor)
    input$labelcolor2
  })
  labelcolorclick <- reactive({
    validcolor <- is_valid_color(trimws(input$labelcolorclick))
    shinyFeedback::feedbackWarning("labelcolorclick", !validcolor, "Please select a valid color or HEX code")
    req(validcolor)
    input$labelcolorclick
  })
  outlinecolor <- reactive({
    validcolor <- is_valid_color(trimws(input$outlinecolor))
    shinyFeedback::feedbackWarning("outlinecolor", !validcolor, "Please select a valid color or HEX code")
    req(validcolor)
    input$outlinecolor
  })
  thresholdcolor <- reactive({
    validcolor <- is_valid_color(trimws(input$thresholdcolor))
    shinyFeedback::feedbackWarning("thresholdcolor", !validcolor, "Please select a valid color or HEX code")
    req(validcolor)
    input$thresholdcolor
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
  
  # Make reactive lists to store specific genes to label
  # These are simply lists that contain gene names. The dataframes that store information about these genes are defined below. 
  genes_to_label_1 <- reactiveVal(c())
  genes_to_label_2 <- reactiveVal(c())
  genes_to_label_click <- reactiveVal(c())
  
  # Reset button for clicked labels
  observeEvent(input$reset_clicked_labels, {
    genes_to_label_click(c())
  })
  
 # Observe plotly click events in which the user clicks on points to add the gene to the clicked gene list
  # Priority must be set to "event" or else Plotly will ignore multiple clicks in the same spot
  observeEvent(event_data("plotly_click", priority = "event"), {
    click_data <- event_data("plotly_click")
    req(click_data)
    
    # Get clicked gene info
    clicked_gene <- data() %>%
      filter(avg_log2FC == click_data$x & neg_log10_pval == click_data$y) %>%
      select(gene) %>%
      pull()

    # Append clicked gene information 
    # If the click successfully pulled a gene:
    if(length(clicked_gene) > 0) { 
      # If the gene is not in genes_to_label_click OR in other gene input lists,
      # Add the gene to genes_to_label_click
      if(!(clicked_gene %in% genes_to_label_click() ||
           clicked_gene %in% genes_to_label_1() ||
           clicked_gene %in% genes_to_label_2())) {
        updated_data <- c(genes_to_label_click(), clicked_gene)
        genes_to_label_click(updated_data)

        # If the gene is already in genes_to_label_click AND NOT in the other lists,
        # remove it: 
      } else if (clicked_gene %in% genes_to_label_click() &&
                 !(clicked_gene %in% genes_to_label_1()) &&
                 !(clicked_gene %in% genes_to_label_2())) {
        updated_data_2 <- setdiff(genes_to_label_click(), clicked_gene)
        genes_to_label_click(updated_data_2)
      }
    }
  }, ignoreNULL = FALSE)
  
   # Process gene input 1
  observeEvent(input$gene_input_1, {
    # Split input by comma and remove empty strings
    input_genes <- unlist(strsplit(input$gene_input_1, "[,\n]+"))
    input_genes <- trimws(input_genes)
    input_genes <- input_genes[input_genes != ""]
    # Update genes_to_label dynamically
    genes_to_label_1(input_genes)
    # If any gene is in genes_to_label_click(), remove it 
    genes_to_label_click(setdiff(genes_to_label_click(), genes_to_label_1()))
  })
  
  # Process gene input 2
  observeEvent(input$gene_input_2, {
    # Split input by comma and remove empty strings
    input_genes <- unlist(strsplit(input$gene_input_2, "[,\n]+"))
    input_genes <- trimws(input_genes)
    input_genes <- input_genes[input_genes != ""]
    # Update genes_to_label dynamically 
    genes_to_label_2(input_genes)
    # If any gene is in genes_to_label_click(), remove it 
    genes_to_label_click(setdiff(genes_to_label_click(), genes_to_label_2()))
  })
  
  # Create reactive dataframe for data, using input data and process_for_volcano function
  de <- reactive({
    de <- NULL
    if(!is.null(data())) {
      de <- process_for_volcano(de = data(),
                                log_fc_cutoff = logfc(),
                                p_val_cutoff = pval(),
                                genes_to_label_1 = genes_to_label_1(),
                                genes_to_label_2 = genes_to_label_2(),
                                genes_to_label_click = genes_to_label_click())
    }
    return(de)
  })
  
  # Initialize annotation dataframes as a list reactive values
  # all_annotations_positions = a dataframe that contains a list of all the genes that are currently annotated with their positions
  # annotations_1 = a dataframe that contains all information for the genes annotated in layer 1
  # annotations_2 = a dataframe that contains all information for the genes annotated in layer 2
  # annotations_click = a dataframe that contains all information for the genes annotated by clicking 
  dataframes <- reactiveValues(
    all_annotations_positions = reactiveVal(data.frame(gene = character(), ax = numeric(), ay = numeric())),
    annotations_1 = reactiveVal(data.frame(gene = character(), 
                                           avg_log2FC = numeric(), 
                                           p_val_adj = numeric(), 
                                           neg_log10_pval = numeric(), 
                                           name1 = character(), 
                                           ax = numeric(), 
                                           ay = numeric())),
    annotations_2 = reactiveVal(data.frame(gene = character(), 
                                           avg_log2FC = numeric(), 
                                           p_val_adj = numeric(), 
                                           neg_log10_pval = numeric(), 
                                           name2 = character(), 
                                           ax = numeric(), 
                                           ay = numeric())),
    annotations_click = reactiveVal(data.frame(gene = character(), 
                                           avg_log2FC = numeric(), 
                                           p_val_adj = numeric(), 
                                           neg_log10_pval = numeric(), 
                                           nameclick = character(), 
                                           ax = numeric(), 
                                           ay = numeric())),
  )
    
  # Create function to update annotations_1
  update_annotations_1 <- reactive({
    # Get genes to label in layer 1 from the de dataframe
    annot_de_1 <- de()[de()$name1 != "", ]
    # Initiate annotations dataframe with default positions
    if(nrow(annot_de_1) > 0) {
       annot_de_1$ax <- 20
       annot_de_1$ay <- -20
       # If any of these genes are already annotated and stored in all_annotations_positions, update annot_de_1 with gene positions
       aap <- dataframes$all_annotations_positions
       if(nrow(aap) > 0) {
         merged1 <- merge(annot_de_1, aap, by.x = "name1", by.y = "gene", all.x = TRUE)
         merged1$ax <- ifelse(is.na(merged1$ax.y), merged1$ax.x, merged1$ax.y)
         merged1$ay <- ifelse(is.na(merged1$ay.y), merged1$ay.x, merged1$ay.y)
         annot_de_1 <- merged1
       }
       # Finish processing annot_de_1 dataframe by sorting and selecting columns of interest 
       annot_de_1 <- annot_de_1 %>% arrange(desc(neg_log10_pval))
       annot_de_1 <- annot_de_1 %>% select(gene, avg_log2FC, p_val_adj, neg_log10_pval, name1, ax, ay)
    }
    # Return updated dataframe
    dataframes$annotations_1 <- annot_de_1
  })
  
  # Create function to update annotations_2
  update_annotations_2 <- reactive({
    # Get genes to label in layer 2 from the de dataframe
    annot_de_2 <- de()[de()$name2 != "", ]
    if(nrow(annot_de_2) > 0) {
      annot_de_2$ax <- 20
      annot_de_2$ay <- -20
      # If any genes are already in all_annotations_positions, update annot_de_2 with gene positions
      aap <- dataframes$all_annotations_positions
      if(nrow(aap) > 0) {
        merged2 <- merge(annot_de_2, aap, by.x = "name2", by.y = "gene", all.x = TRUE)
        merged2$ax <- ifelse(is.na(merged2$ax.y), merged2$ax.x, merged2$ax.y)
        merged2$ay <- ifelse(is.na(merged2$ay.y), merged2$ay.x, merged2$ay.y)
        annot_de_2 <- merged2
      }
      # Finish processing annot_de_2 dataframe
      annot_de_2 <- annot_de_2 %>% arrange(desc(neg_log10_pval))
      annot_de_2 <- annot_de_2 %>% select(gene, avg_log2FC, p_val_adj, neg_log10_pval, name2, ax, ay)
    }
    # Return updated dataframe
    dataframes$annotations_2 <- annot_de_2
  })
  
  # Create function to update annotations_click
  update_annotations_click <- reactive({
    # Get genes to label in "click" layer from the de dataframe
    annot_de_click <- de()[de()$nameclick != "", ]
    if(nrow(annot_de_click) > 0) {
      annot_de_click$ax <- 20
      annot_de_click$ay <- -20
    # If there are already genes in all_annotations_positions, update annot_de_click with gene positions
      aap <- dataframes$all_annotations_positions
      if(nrow(aap) > 0) {
        mergedclick <- merge(annot_de_click, aap, by.x = "nameclick", by.y = "gene", all.x = TRUE)
        mergedclick$ax <- ifelse(is.na(mergedclick$ax.y), mergedclick$ax.x, mergedclick$ax.y)
        mergedclick$ay <- ifelse(is.na(mergedclick$ay.y), mergedclick$ay.x, mergedclick$ay.y)
        annot_de_click <- mergedclick
      }
      # Finish processing annot_de_click dataframe
      annot_de_click <- annot_de_click %>% arrange(desc(neg_log10_pval))
      annot_de_click <- annot_de_click %>% select(gene, avg_log2FC, p_val_adj, neg_log10_pval, nameclick, ax, ay)
    }
    # Return updated dataframe
    dataframes$annotations_click <- annot_de_click
  })
  
  # Create function to update all_annotations_positions from the data in the annotations_1, annotations_2, and annotations_click dataframes
  update_all_annotations_positions <- reactive({
    aap <- data.frame(gene = c(), ax = numeric(), ay = numeric())
    if(nrow(dataframes$annotations_1) > 0 || nrow(dataframes$annotations_2) > 0 || nrow(dataframes$annotations_click > 0)) {
      aap <- bind_rows(dataframes$annotations_1, dataframes$annotations_2, dataframes$annotations_click)
      aap <- aap %>% select(gene, ax, ay)
    }
    dataframes$all_annotations_positions <- aap
  })
  
  # Generate the plot using the make_volcano_plotly function
  plot <- reactive({
    p <- NULL
    if(!is.null(de())) {
      p <- make_volcano_plotly(
        de(),
        dataframes$annotations_1,
        dataframes$annotations_2,
        dataframes$annotations_click,
        graph_title = title(),
        log_fc_cutoff = logfc(), 
        p_val_cutoff = pval(),
        upcolor = upcolor(), 
        downcolor = downcolor(), 
        midcolor = midcolor(),
        labelcolor1 = labelcolor1(),
        labelsize1 = labelsize1(),
        labelcolor2 = labelcolor2(),
        labelsize2 = labelsize2(),
        labelcolorclick = labelcolorclick(),
        labelsizeclick = labelsizeclick(),
        pointsize = pointsize(),
        visible_cutoffs = show_cutoffs(),
        show_outlines = show_outlines(),
        outline_color = outlinecolor(),
        height = height(),
        width = width(),
        space = space(),
        threshold_color = thresholdcolor(),
        download.filetype = filetype())
    }
    return(p)
  })
  
  
  # Render the plot
  output$plot <- renderPlotly(
    {
      update_annotations_1()
      update_annotations_2()
      update_annotations_click()
      update_all_annotations_positions()
      plot()
    
    })
  
  # Observe events
  # This is where we check to see if the user has clicked and dragged a label to reposition it
  observeEvent(event_data("plotly_relayout"), {
    relayout_data <- event_data("plotly_relayout")
    req(relayout_data)
    
    # Check if names in relayout_data start with "annotations"
    annotation_names <- names(relayout_data)
    annotation_indices <- grep("^annotations\\[[0-9]+\\]", annotation_names)
    
    if(length(annotation_indices) > 0) {
      for(index in annotation_indices) {
        name <- annotation_names[index]
        # Extract index number
        i <- as.numeric(sub("annotations\\[([0-9]+)\\]\\..*", "\\1", name))
        
        # Extract ax and ay values
        if(grepl("\\.ax$", name)) {
          ax_value <- relayout_data[[name]]
        } else if (grepl("\\.ay$", name)) {
          ay_value <- relayout_data[[name]]
        }
      }
      # Update all_annotations_positions with new position info
      aap <- dataframes$all_annotations_positions
      aap$ax[i + 1] <- ax_value
      aap$ay[i + 1] <- ay_value
      dataframes$all_annotations_positions <- aap
    } 
  })
  
  
} # End of server function
  
    

# Run the application 
shinyApp(ui = ui, server = server)
