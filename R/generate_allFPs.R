#' Generate Feature Plots for given genes across different clusters
#'
#' This function generates feature plots for each gene provided in the gene_list parameter across
#' different clusters in a Seurat object. It arranges these plots and saves them as a PDF file.
#'
#' @param SeuFile A Seurat object that contains the gene expression data.
#' @param gene_list A list of genes for which the feature plots need to be generated.
#'
#' @return Generates a PDF file containing the feature plots for the provided genes.
#'
#' @importFrom Seurat DefaultAssay FeaturePlot
#' @importFrom ggpubr ggarrange
#' @importFrom grDevices pdf dev.off
#'
#' @export


generate_AllFPs <- function(SeuFile, gene_list){

  Filename = as.character(substitute(SeuFile))

  DefaultAssay(SeuFile) <- "RNA"

  # Cleaned gene_list
  gene_list_2 = list()
  for(d in names(gene_list)){gene_list_2[[d]] = subset(gene_list[[d]], gene_list[[d]] %in% row.names(SeuFile))}

  FPList = list()

  for(d in names(gene_list_2)){
    Genes = gene_list_2[[d]]

    if(length(Genes) > 1){
      FPSinglePage = list()
      FPSinglePage[[1]] = FeaturePlot(SeuFile, Genes[1], reduction="umap") + labs(title=paste(d, Genes[1]))
      for(p in seq(2, length(Genes), 1)){
        FPSinglePage[[p]] = FeaturePlot(SeuFile, Genes[p], reduction="umap")
      }
      FPList[[d]] = ggarrange(plotlist = FPSinglePage, ncol=5, nrow=4)
    } else{
      FPSinglePage = list()
      FPSinglePage[[1]] = FeaturePlot(SeuFile, Genes[1], reduction="umap") + labs(title=paste(d, Genes[1]))
      FPList[[d]] = ggarrange(plotlist = FPSinglePage, ncol=5, nrow=4)
    }
  }

  if (!file.exists("output")) {
    dir.create("output")
  }

  pdf(paste("output/", Filename, "_AllFPs.pdf", sep=""),
      width = 25, height = 22)
  print(FPList)
  dev.off()

}
