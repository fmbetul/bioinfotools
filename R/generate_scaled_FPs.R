



generate_scaled_FPs <- function(seurat_objects, genes, filename, seperate = T) {

  # Find the maximum value for each gene across all Seurat objects
  max_values <- sapply(genes, function(gene) {
    max(sapply(seurat_objects, function(seurat_obj) {

      # Set the default assay for the Seurat object to 'RNA'
      DefaultAssay(seurat_obj) <- 'RNA'

      # Check if the gene is present in the RNA assay
      if (gene %in% rownames(seurat_obj@assays$RNA@counts)) {

        # Fetch the gene expression data
        max(Seurat::FetchData(seurat_obj, vars = gene), na.rm = TRUE)
      } else {
        NA
      }
    }))
  }, USE.NAMES = TRUE)



  # Loop through the Seurat objects
  for (seurat_obj_name in names(seurat_objects)) {

    # Get the current Seurat object
    seurat_obj <- seurat_objects[[seurat_obj_name]]

    # Set the default assay for the Seurat object to 'RNA'
    DefaultAssay(seurat_obj) <- 'RNA'

    if(length(unique(genes)) == 1){setncol = 1; setnrow =1
    }else if(length(unique(genes)) > 1 & length(unique(genes)) <= 2){setncol = 2; setnrow =1
    }else if(length(unique(genes)) > 2 & length(unique(genes)) <= 4){setncol = 2; setnrow =2
    }else if(length(unique(genes)) > 4 & length(unique(genes)) <= 9){setncol = 3; setnrow =3
    }else if(length(unique(genes)) > 9 & length(unique(genes)) <= 16){setncol = 4; setnrow =4
    }else if(length(unique(genes)) > 16){setncol = 4; setnrow =ceiling(length(unique(genes))/4)
    }

    # Create a combined feature plot for all genes in the Seurat object
    p <- FeaturePlot(seurat_obj, ncol= setncol, features = genes, combine = TRUE)

    # Loop through the genes and add a color scale to each feature plot
    for (i in seq_along(genes)) {
      gene <- genes[[i]]
      p[[i]] <- p[[i]] + scale_color_gradientn(colours = c("lightgrey", "navy"), limits = c(0, max_values[[gene]]))
    }

    # Save the feature plots as a PDF file
    if (!file.exists("output")) {
      dir.create("output")
    }

    if(seperate == T){
      pdf_filename <- paste0("output/", filename, "_", seurat_obj_name, "_FeaturePlots.pdf")
      pdf(file = pdf_filename, width = 3.75*setncol, height = 3*setnrow)
      print(p)
      dev.off()
    }else if(seperate == F){
      if(seurat_obj_name == names(seurat_objects)[[1]]){
        CombinePlots = list()
        CombinePlots[[seurat_obj_name]] = p + plot_annotation(title = seurat_obj_name)
      }else{
        CombinePlots[[seurat_obj_name]] = p + plot_annotation(title = seurat_obj_name)
      }

    }
  }
  if(seperate == F){
    pdf_filename <- paste0("output/", filename, "_FeaturePlots.pdf")
    pdf(file = pdf_filename, width = 3.75*setncol, height = 3*setnrow)
    print(CombinePlots)
    dev.off()
  }
}
