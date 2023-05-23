#' Plot specific cells from an integrated Seurat object
#'
#' This function generates a UMAP plot for specific cells in a given sample that
#' belong to a certain cluster from an integrated Seurat object.
#'
#' @param integrated A Seurat object. This is an integrated or merged Seurat object
#'                  that includes data from multiple samples.
#' @param ClusterName A character string specifying the name of the cell type or
#'                   cluster whose cells are to be visualized.
#' @param sample A Seurat object that contains scRNA-seq data for a specific
#'               sample.
#'
#' @return This function does not return a value, but produces a UMAP plot for
#'         the selected cells from the integrated Seurat object.
#'
#' @seealso \code{\link[Seurat]{DimPlot}}
#'
#' @importFrom Seurat DimPlot
#'
#' @export


plot_cells_from_integrated <- function(integrated, ClusterName, sample){

  sample_name <- as.character(substitute(sample))

  # Get the cell barcodes for the selected cluster in integrated
  cell.barcodes <- rownames(subset(integrated@meta.data, integrated@meta.data$celltype == ClusterName))

  # Select the cell barcodes that start with sample name
  cell.barcodes <- grep(paste0("^", sample_name, "_"), cell.barcodes, value = TRUE)

  # Remove the sample prefix from the selected cell barcodes
  cell.barcodes <- sub(paste0("^", sample_name, "_"), "", cell.barcodes)

  # Visualize the selected cluster cells from integrated
  DefaultAssay(sample) <- "RNA"

  if(is.null(sample@meta.data$celltype) == FALSE) {

    DimPlot(sample, group.by = 'celltype', reduction = "umap", label = TRUE, cells.highlight = cell.barcodes)

  }else{

    DimPlot(sample, label = TRUE, cells.highlight = cell.barcodes)

  }

}
