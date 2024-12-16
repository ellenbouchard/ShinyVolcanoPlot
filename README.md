# ShinyVolcanoPlot
My own Shiny app for making custom volcano plots. 


My goal is to create an interactive web app for easy volcano plot customization. 
Currently, I have created two apps:
1. A basic Shiny app that allows the user to upload a .csv file of their differential expression results, and lets them adjust thresholds, colors, and perform some label customization via ggplot
2. A Shiny app that renders the plot in Plotly, allowing for more fine-tuned customizaitons, including drag-and-drop label positioning.

Currently, I am continuing to work on the Plotly app to add the following functionality:

- Cutsomize appearance of specific gene markers (point color, size, outline)
- Download function
