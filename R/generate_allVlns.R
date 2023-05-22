#' Generate Violin Plots for given genes across different clusters
#'
#' This function generates violin plots for each gene provided in the gene_list parameter across
#' different clusters in a Seurat object. It arranges these plots and saves them as a PDF file.
#'
#' @param SeuFile A Seurat object that contains the gene expression data.
#' @param gene_list A list of genes for which the violin plots need to be generated.
#'
#' @return Generates a PDF file containing the violin plots for the provided genes.
#'
#' @importFrom Seurat DefaultAssay VlnPlot
#' @importFrom ggpubr ggarrange
#' @importFrom grDevices pdf dev.off
#' @importFrom ggplot2 coord_flip theme element_blank element_text labs
#'
#' @export


generate_allVlns <- function(SeuFile, gene_list){

  Filename = as.character(substitute(SeuFile))

  DefaultAssay(SeuFile) <- "RNA"


  # Cleaned gene_list
  gene_list_2 = list()
  for(d in names(gene_list)){gene_list_2[[d]] = subset(gene_list[[d]], gene_list[[d]] %in% row.names(SeuFile))}

  AllVlns = list()

  for(d in names(gene_list_2)){
    Genes = gene_list_2[[d]]
    StorePlots = list()

    if(length(Genes) > 1){

      for(x in Genes[1]){
        plotA <- VlnPlot(SeuFile, features = x, pt.size = 0, same.y.lims = F,)
        plotA <- plotA + coord_flip()+ theme(axis.ticks.x= element_blank(), axis.title.y=element_blank(), axis.text.x=element_blank(), axis.title.x=element_blank(),
                                             axis.ticks.y = element_blank(), legend.position = "none", plot.title = element_text(size=12))+ labs(title = d, subtitle = Genes[1])
        StorePlots[[x]] = plotA
      }

      for(x in Genes[2:length(Genes)]){
        plotB <- VlnPlot(SeuFile, features = x, pt.size = 0, same.y.lims = F,)
        plotB <- plotB + coord_flip()+ theme(axis.ticks.x= element_blank(), axis.title.y=element_blank(), axis.text.x=element_blank(), axis.title.x=element_blank(),
                                             axis.ticks.y = element_blank(), legend.position = "none", axis.text.y = element_blank(), plot.title = element_text(size=12))
        StorePlots[[x]] = plotB
      }
      AllVlns[[d]] <- ggarrange(plotlist = StorePlots, widths=c(1.4, rep(1, length(Genes)-1)), ncol = 20,  nrow = 1)

    } else {

      for(x in Genes[1]){
        plotA <- VlnPlot(SeuFile, features = x, pt.size = 0, same.y.lims = F,)
        plotA <- plotA + coord_flip()+ theme(axis.ticks.x= element_blank(), axis.title.y=element_blank(), axis.text.x=element_blank(), axis.title.x=element_blank(),
                                             axis.ticks.y = element_blank(), legend.position = "none", plot.title = element_text(size=12))+ labs(title = d, subtitle = Genes[1])
        StorePlots[[x]] = plotA
      }
      AllVlns[[d]] <- ggarrange(plotlist = StorePlots, widths=c(1.4, rep(1, length(Genes)-1)), ncol = 20,  nrow = 1)
    }
  }


  if (!file.exists("output")) {
    dir.create("output")
  }

  pdf(paste("output/", Filename, "_AllVlns.pdf", sep=""),
      width=40, height=length(unique(SeuFile@active.ident)))
  print(AllVlns)
  dev.off()
}
