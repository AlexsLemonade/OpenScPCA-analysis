# barchart with or without faceting
# each bar is a stacked barchart using the fill_color
# faceting is only done if a facet_variable is provided
stacked_barchart <- function(
    df, 
    fill_column,
    celltype_colors, # named vector where names match the values in fill_column
    facet_variable = NULL
){
  
  barchart <- ggplot(df) + 
    aes(
      x = library_id, 
      y = percent_cells_annotation, 
      fill = !!sym(fill_column)
    ) +
    geom_col() + 
    scale_y_continuous(expand = c(0,0)) +
    scale_fill_manual(values = celltype_colors) +
    theme(axis.text.x = element_blank()) +
    labs(
      fill= "cell type"
    )
  
  if(!is.null(facet_variable)){
    barchart <- barchart +
      facet_wrap(vars(!!sym(facet_variable)), scales ="free_x")
  }
  
  return(barchart)
}

# sina plot looking at immune percentage on the y-axis 
sina_plot <- function(df, plot_column){
  
  ggplot(df, aes(x = !!sym(plot_column), y = percent_immune)) +
    ggforce::geom_sina() +
    stat_summary(fun.y=median, geom="crossbar" , color = "red", linewidth = 0.2) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          plot.margin = margin(10,10,10,10)) +
    labs(
      x = "", 
      y = "Percent of cells annotated as Immune",
      title = plot_column
    )
  
}
