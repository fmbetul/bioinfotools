#' Plot cells of a specified sample and cluster type
#'
#' This function plots cells of a specified sample and cluster type using a UMAP plot.
#'
#' @param sample A Seurat object that contains the sample data.
#' @param ClusterName A string that specifies the cell type or cluster name to be plotted.
#' @param integrated An integrated Seurat object that contains the merged data from multiple samples.
#'
#' @return The function does not return anything. It prints a UMAP plot to the console.
#'
#' @seealso \code{\link[Seurat]{DimPlot}}
#'
#' @importFrom Seurat DimPlot
#'
#' @export


plot_cells_from_sample <- function(sample, ClusterName, integrated){

  sample_name <- as.character(substitute(sample))

  # Get the cell barcodes
  cell.barcodes <- rownames(subset(sample@meta.data, sample@meta.data$celltype == ClusterName))

  # Add the sample prefix
  cell.barcodes <- paste0(sample_name, "_", cell.barcodes)

  # Visualize the Cluster cells from sample Seurat object
  print(DimPlot(integrated, group.by = 'celltype', reduction = "umap", label = TRUE, cells.highlight = cell.barcodes))

}

