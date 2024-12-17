#
# This is my pet project: a Shiny app that renders custom, dynamic volcano plots
# Date Created: December 4 2024
#

library(shiny)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(plotly)
library(shinyjs)

### FUNCTIONS

## Function 1: process data for plotting
# (This used to be in the same function as make_volcano, but for the purposes
# of dynamic annotations needs to be separate)

process_for_volcano = function(de,
                             log_fc_cutoff = 0.05,
                             p_val_cutoff = 0.01,
                             genes_to_label_1 = c(),
                             genes_to_label_2 = c(),
                             genes_to_label_click = c()
                             ) {
  
  # Set cutoffs
  up_cutoff <- log_fc_cutoff
  down_cutoff <- -1 * log_fc_cutoff
  
  # For adjusted p values of zero, change neg_log10_pval so that it is not infinite
  de$neg_log10_pval[de$p_val_adj == 0] <- -log10(.Machine$double.xmin)
  
  # Label genes as significantly "UP" or "DOWN" based on p val and logFC thresholds
  de$reg = ""
  de$reg[de$p_val_adj < p_val_cutoff & de$avg_log2FC > up_cutoff & de$avg_log2FC > 0] <- "UP"
  de$reg[de$p_val_adj < p_val_cutoff & de$avg_log2FC < down_cutoff & de$avg_log2FC < 0] <- "DOWN"
  
  # Add first layer of labels
  de$name1 = ""
  if(length(genes_to_label_1) > 0) {de <- de %>% mutate(name1 = ifelse(gene  %in% genes_to_label_1 & reg != "",gene, ""))}  
  
  # Add second layer of labels
  de$name2 = ""
  if(length(genes_to_label_2) > 0) {de <- de %>% mutate(name2 = ifelse(gene  %in% genes_to_label_2 & reg != "",gene, ""))}  
  
  # Add third ("clicked") layer of labels, only if genes are not in first two layers of labels
  de$nameclick = ""
  if(length(genes_to_label_click) > 0) {
    de <- de %>% mutate(nameclick = ifelse(gene %in% genes_to_label_click & 
                                             reg != "" &
                                             !(gene %in% genes_to_label_1 | gene %in% genes_to_label_2), 
                                           gene, ""))
    }
  
  return(de)
}

# Function make_volcano_ploty: takes the output of process_for_volcano and returns volcano plot made with Plotly
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
  
  up_cutoff <- log_fc_cutoff
  down_cutoff <- -1 * log_fc_cutoff

  # Define colors
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

  
  # Add black out lines if show_outlines is TRUE
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
    
  } # Add markers
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
                          dragmode = FALSE)
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
                 selectInput("filetype", "Select File Type", choices = c("png","svg"), selected = "png"),
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


# SERVER FUNCTION
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
  
  # Define input variables as reactives
  title <- debounce(reactive(input$title), millis = 300)
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
  genes_to_label_1 <- reactiveVal(c())
  genes_to_label_2 <- reactiveVal(c())
  genes_to_label_click <- reactiveVal(c())
  
  # Make reactive dataframe to store just position information (ax and ay) of all annotation layers
  # This is necessary for retaining position of labels that are repositioned by the user 
  all_annotations_positions <- reactiveVal(data.frame(gene = character(), ax = numeric(), ay = numeric()))
  
 # Observe plotly click events
  # priority must be set to "event" or else Plotly will ignore multiple clicks in the same spot
  observeEvent(event_data("plotly_click", priority = "event"), {
    click_data <- event_data("plotly_click")
    req(click_data)
    
    # Get clicked gene info
    clicked_gene <- data() %>%
      filter(avg_log2FC == click_data$x & neg_log10_pval == click_data$y) %>%
      select(gene) %>%
      pull()

    # Append clicked gene information 
    if(length(clicked_gene) > 0) { # If the click successfully pulled a gene:
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
  
  # Render the plot
  output$plot <- renderPlotly(
    {
      # First, process data
      # This has to be done every time to account for multiple annotation layers
    de <- process_for_volcano(de = data(),
                                log_fc_cutoff = logfc(),
                                p_val_cutoff = pval(),
                                genes_to_label_1 = genes_to_label_1(),
                                genes_to_label_2 = genes_to_label_2(),
                                genes_to_label_click = genes_to_label_click())
    
    # Generate dataframes for lists of labeled genes
    annot_de_1 <- de[de$name1 != "", ]
    annot_de_2 <- de[de$name2 != "", ]
    annot_de_click <- de[de$nameclick != "", ]
    
    # Start by populating annotation dataframes with default label positions (ax  = 20 and ay = -20)
    # Then, merge each annotation dataframe with all_annotations_positison()
    # such that if a label already has a non-default position on the plot, it's retained
    if(nrow(annot_de_1) > 0) {
      annot_de_1$ax <- 20
      annot_de_1$ay <- -20
      merged1 <- merge(annot_de_1, all_annotations_positions(), by.x = "name1", by.y = "gene", all.x = TRUE)
      merged1$ax <- ifelse(is.na(merged1$ax.y), merged1$ax.x, merged1$ax.y)
      merged1$ay <- ifelse(is.na(merged1$ay.y), merged1$ay.x, merged1$ay.y)
      annot_de_1 <- merged1 %>% arrange(desc(neg_log10_pval))
    }
    if(nrow(annot_de_2) > 0) {
      annot_de_2$ax <- 20
      annot_de_2$ay <- -20
      merged2 <- merge(annot_de_2, all_annotations_positions(), by.x = "name2", by.y = "gene", all.x = TRUE)
      merged2$ax <- ifelse(is.na(merged2$ax.y), merged2$ax.x, merged2$ax.y)
      merged2$ay <- ifelse(is.na(merged2$ay.y), merged2$ay.x, merged2$ay.y)
      annot_de_2 <- merged2 %>% arrange(desc(neg_log10_pval))
    }
    if(nrow(annot_de_click) > 0) {
      annot_de_click$ax <- 20
      annot_de_click$ay <- -20
      mergedclick <- merge(annot_de_click, all_annotations_positions(), by.x = "nameclick", by.y = "gene", all.x = TRUE)
      mergedclick$ax <- ifelse(is.na(mergedclick$ax.y), mergedclick$ax.x, mergedclick$ax.y)
      mergedclick$ay <- ifelse(is.na(mergedclick$ay.y), mergedclick$ay.x, mergedclick$ay.y)
      annot_de_click <- mergedclick %>% arrange(desc(neg_log10_pval))
    }
    
    # Update all_annotations_positions with newly added annotations
    if(nrow(annot_de_1) > 0 || nrow(annot_de_2) > 0 || nrow(annot_de_click) > 0) {
      aap <- bind_rows(annot_de_1, annot_de_2, annot_de_click)
      aap <- aap %>% select(gene, ax, ay)
      all_annotations_positions(aap)
    }
    
    # Finally, make the plot
    p <- make_volcano_plotly(
                 de,
                 annot_de_1,
                 annot_de_2,
                 annot_de_click,
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
    return(p)
    
    })
  
  # Observe events
  # This is where we check to see if the user has clicked and dragged a label to reposition
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
          all_annotations_positions()
        } else if (grepl("\\.ay$", name)) {
          ay_value <- relayout_data[[name]]
        }
      }
      # Update all_annotations_positions with new position info
      aap <- all_annotations_positions()
      aap$ax[i + 1] <- ax_value
      aap$ay[i + 1] <- ay_value
      all_annotations_positions(aap)
    } 
  })
  
  # Render the figure for download
    output$download <- downloadHandler(
      filename = function() {
        #paste("plot", Sys.Date(), ".", input$filetype, sep = "")
        paste("plot", Sys.Date(), ".html", sep = "")
      },
      content = function(file) {
  ## CODE TO SAVE PLOT
        
      })
  
} # End of server function
  
    

# Run the application 
shinyApp(ui = ui, server = server)
