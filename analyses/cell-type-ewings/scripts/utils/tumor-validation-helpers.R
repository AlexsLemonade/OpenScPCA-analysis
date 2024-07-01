# These are functions used to prepare data frames for plotting when doing tumor cell validataion 

# Data frame prep --------------------------------------------------------------

create_classification_df <- function(
    sce,
    marker_gene_results_file, # path to marker gene classifications
    copykat_predictions_file, # path to copykat predictions
    infercnv_predictions_file, # path to infercnv predictions
    infercnv_metadata_file, # path to full infercnv metadata
    geneset_scores_file # path to gene set scores results
){
 
  # read in marker gene results 
  marker_gene_results_df <- readr::read_tsv(marker_gene_results_file, show_col_types = FALSE)
  
  # copykat output
  copykat_results_df <- readr::read_tsv(copykat_predictions_file, show_col_types = FALSE) |> 
    # only select no_ref classification and rename to be consistent with naming of other classifications
    dplyr::select(barcodes, copykat_classification = no_ref, mean_cnv_detection) |> 
    dplyr::distinct() |> 
    dplyr::mutate(
      copykat_classification = dplyr::case_when(
        copykat_classification == "aneuploid" ~ "Tumor",
        copykat_classification == "diploid" ~ "Normal",
        .default = copykat_classification
      )
    )
  
  # infercnv output 
  infercnv_results_df <- readr::read_tsv(infercnv_predictions_file, show_col_types = FALSE)
  # use read.delim to ensure row names do not get read in as a column 
  infercnv_metadata_df <- read.delim(infercnv_metadata_file) |> 
    tibble::rownames_to_column("barcodes")
  
  # remove MT and GL chr, comments in docs recommend that these are not useful and can be inaccurate 
  chr_columns_to_keep <- which(!stringr::str_detect(names(infercnv_metadata_df), "chrMT|chrGL"))
  infercnv_metadata_df <- infercnv_metadata_df[, chr_columns_to_keep]
  
  # gene set scores 
  geneset_scores_df <- readr::read_tsv(geneset_scores_file, show_col_types = FALSE)
  
  # pull out the UMAP coordinates and make a data frame to use for plotting
  classification_df <- sce |> 
    scuttle::makePerCellDF(use.dimred = "UMAP") |> 
    # replace UMAP.1 with UMAP1
    dplyr::rename_with(
      \(x) stringr::str_replace(x, "^UMAP\\.", "UMAP")
    ) |>
    # get rid of excess columns
    dplyr::select(barcodes, UMAP1, UMAP2) |> 
    # add in marker gene classifications
    dplyr::left_join(marker_gene_results_df, by = "barcodes") |> 
    # add in copykat results 
    dplyr::left_join(copykat_results_df, by = "barcodes") |> 
    # add in infercnv results and metadata
    dplyr::left_join(infercnv_results_df, by = "barcodes") |> 
    dplyr::left_join(infercnv_metadata_df, by = "barcodes") |> 
    # gene set scores 
    dplyr::left_join(geneset_scores_df, by = "barcodes")
  
  return(classification_df)
   
}


# create a table with one row per marker gene per cell 
# include sum of the scaled gene expression for all marker genes in a cell 
create_marker_gene_df <- function(
    sce, # sce object to grab gene expression from 
    classification_df, # data frame output from create_classification_df()
    marker_genes_file # path to file containing list of marker genes used
){
  # read in marker genes table 
  marker_genes_df <- readr::read_tsv(marker_genes_file, show_col_types = FALSE) |> 
    # account for genes being from multiple sources
    dplyr::select(cell_type, ensembl_gene_id, gene_symbol) |> 
    dplyr::distinct()
  
  # get list of marker genes to use 
  marker_genes <- marker_genes_df |> 
    dplyr::filter(cell_type == "tumor") |> 
    dplyr::pull(ensembl_gene_id)
  
  # get the gene expression counts for all marker genes
  marker_gene_exp <- logcounts(sce[marker_genes, ]) |>
    as.matrix() |>
    t() |>
    as.data.frame() |>
    tibble::rownames_to_column("barcodes")
  
  # calculate sum of the scaled gene expression values for plotting 
  plot_markers_df <- classification_df |> 
    # add in marker gene expression to dataframe 
    dplyr::left_join(marker_gene_exp, by = "barcodes") |> 
    # combine all genes into a single column for easy faceting
    tidyr::pivot_longer(
      cols = starts_with("ENSG"),
      names_to = "ensembl_gene_id",
      values_to = "gene_expression"
    ) |> 
    # join with marker gene df to get gene symbols for plotting 
    dplyr::left_join(marker_genes_df, by = c("ensembl_gene_id")) |> 
    # keep only gene symbols, barcodes, and classification columns 
    dplyr::select(barcodes, gene_symbol, gene_expression, ends_with("classification")) |> 
    unique() |> 
    dplyr::group_by(gene_symbol) |>
    # get z-scores for each gene
    dplyr::mutate(transformed_gene_expression = scale(gene_expression)[, 1]) |> 
    dplyr::ungroup() |> 
    dplyr::group_by(barcodes) |> 
    # add all the scores together 
    dplyr::mutate(sum_raw_exp = sum(gene_expression), 
                  sum_transformed_exp = sum(transformed_gene_expression)) |>
    dplyr::ungroup()
  
  return(plot_markers_df)
}

# Heatmaps ---------------------------------------------------------------------

# create a heatmap where rows are marker genes or gene sets and columns are cells 
# this function adds a column annotation since cells are columns 
plot_gene_heatmap <- function(
    df, 
    row_title = "",
    legend_title = "", 
    annotation = NULL
){
  
  # plot heatmap of marker genes 
  heatmap <- ComplexHeatmap::Heatmap(
    df,
    # set the color scale based on min and max values
    col = circlize::colorRamp2(seq(min(df), max(df), length = 2), colors = c("white", "#00274C")),
    border = TRUE,
    ## Row parameters
    cluster_rows = TRUE,
    row_title = row_title,
    row_title_gp = grid::gpar(fontsize = 12),
    row_title_side = "left",
    row_names_side = "left",
    row_dend_side = "right",
    row_names_gp = grid::gpar(fontsize = 10),
    ## Column parameters
    cluster_columns = TRUE,
    show_column_names = FALSE,
    bottom_annotation = annotation,
    heatmap_legend_param = list(
      title = legend_title
    )
  )
  
  return(heatmap)
}

# create a heatmap where rows are cells and columns are chromosomes 
# this function adds a row annotation since cells are rows 
plot_cnv_heatmap <- function(
    df,
    cnv_col_prefix, # prefix used to grab columns with cnv proportion (e.g., "proportion_scaled_loss")
    annotation = NULL,
    legend_title = ""
){
  
  # barcode by chromosome df 
  cnv_df <- df |> 
    dplyr::select(barcodes, starts_with(cnv_col_prefix)) |> 
    unique() |> 
    tibble::column_to_rownames("barcodes")
  
  # just get the name of the chromosome as column name 
  colnames(cnv_df) <- stringr::word(colnames(cnv_df), -1, sep = "_")
  
  # do a little bit of column reordering before converting to a matrix
  cnv_df <- cnv_df |> 
    dplyr::select(paste0("chr", 1:22)) |> 
    as.matrix()
  
  heatmap <- ComplexHeatmap::Heatmap(
    cnv_df,
    col = circlize::colorRamp2(seq(min(cnv_df), max(cnv_df), length = 2), colors = c("white", "#00274C")),
    border = TRUE,
    ## Row parameters
    cluster_rows = TRUE,
    row_title = "",
    row_title_gp = grid::gpar(fontsize = 12),
    row_title_side = "left",
    show_row_names = FALSE,
    row_dend_side = "right",
    row_names_gp = grid::gpar(fontsize = 10),
    ## Column parameters
    cluster_columns = FALSE,
    show_column_names = TRUE,
    right_annotation = annotation,
    heatmap_legend_param = list(
      title = legend_title
    )
  )
  return(heatmap)
}
