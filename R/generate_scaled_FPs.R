#' Generate Scaled Feature Plots
#'
#' This function generates feature plots for a list of genes in multiple Seurat objects,
#' scaling the plots by the maximum value of each gene across all objects.
#' The generated plots are then saved as PDF files in an output directory.
#'
#' @param seurat_objects A named list of Seurat objects.
#' @param genes A character vector of gene names for which to generate feature plots.
#' @param filename The filename (without extension) to use for the output PDF files.
#' @param separate A boolean indicating whether to save a separate PDF for each Seurat object.
#'                Default is TRUE.
#'
#' @return This function does not return a value. It generates PDF files with feature plots.
#'
#' @importFrom Seurat DefaultAssay FeaturePlot FetchData
#' @importFrom ggplot2 scale_color_gradientn element_blank labs ggtitle
#' @importFrom ggpubr ggarrange
#' @importFrom grDevices pdf dev.off
#'
#' @export


generate_scaled_FPs <- function (seurat_objects, genes, filename = as.character(Sys.Date()), separate = F) {

  # Find the maximum value for each gene across all Seurat objects
  max_values <- sapply(genes, function(gene) {
    max(sapply(seurat_objects, function(seurat_obj) {

      # Set the default assay for the Seurat object to 'RNA'
      DefaultAssay(seurat_obj) <- "RNA"

      # Check if the gene is present in the RNA assay
      if (gene %in% rownames(seurat_obj@assays$RNA@counts)) {

        # Fetch the gene expression data
        max(Seurat::FetchData(seurat_obj, vars = gene), na.rm = TRUE)
      } else {
        NA
      }
    }))}, USE.NAMES = TRUE)


  # Loop through the Seurat objects
  for (seurat_obj_name in names(seurat_objects)) {

    # Get the current Seurat object
    seurat_obj <- seurat_objects[[seurat_obj_name]]

    # Set the default assay for the Seurat object to 'RNA'
    DefaultAssay(seurat_obj) <- "RNA"

    # set ncol and nrow
    if (length(unique(genes)) == 1) {
      setncol = 1
      setnrow = 1
    } else if (length(unique(genes)) > 1 & length(unique(genes)) <= 2) {
      setncol = 2
      setnrow = 1
    } else if (length(unique(genes)) > 2 & length(unique(genes)) <= 4) {
      setncol = 2
      setnrow = 2
    } else if (length(unique(genes)) > 4 & length(unique(genes)) <= 9) {
      setncol = 3
      setnrow = 3
    } else if (length(unique(genes)) > 9 & length(unique(genes)) <= 16) {
      setncol = 4
      setnrow = 4
    } else if (length(unique(genes)) > 16) {
      setncol = 4
      setnrow = ceiling(length(unique(genes))/4)
    }


    # Generate a combined feature plot for all genes in the Seurat object
    p <- FeaturePlot(seurat_obj, ncol = setncol, features = genes, combine = TRUE)


    # Loop through the genes and add a color scale to each feature plot
    for (i in seq_along(genes)) {
      gene <- genes[[i]]
      p[[i]] <- p[[i]] + scale_color_gradientn(colours = c("lightgrey", "navy"),
                                               limits = c(0, max_values[[gene]]))

      # Add Seurat object name as title before the name of the first gene
      if (i == 1) {
        p[[i]] <- p[[i]] + ggtitle(paste(seurat_obj_name, genes[1]))
      }
    }

    # Save the feature plots as a PDF file if separate = T
    # Generate CombinePlots object if separate = F
    if (separate == T) {
      pdf_filename <- paste0(filename, "_", seurat_obj_name, "_FPs.pdf")
      pdf(file = pdf_filename, width = 3.75 * setncol, height = 3 * setnrow)
      print(p)
      dev.off()
    }
    else {
      if (seurat_obj_name == names(seurat_objects)[[1]]) {
        CombinePlots = list()
        CombinePlots[[seurat_obj_name]] = p
      }
      else {
        CombinePlots[[seurat_obj_name]] = p
      }
    }
  } # end of for loop for seurat_objects

  # Save the feature plots as a PDF file if separate = F
  if (separate == F) {
    pdf_filename <- paste0(filename, "_FPs.pdf")
    pdf(file = pdf_filename, width = 3.75 * setncol, height = 3 * setnrow)
    print(CombinePlots)
    dev.off()
  }
}
