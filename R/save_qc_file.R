#' Save QC file
#'
#' This function saves UMAP plot and quality control (QC) files, with options to include QC plots and heatmaps, to output folder.
#'
#' @param SeuFile Seurat object that contains the data.
#' @param version A string that provides the version number for the QC file. Default is "v1.0.0".
#' @param qc.plot Logical, if TRUE the function will generate QC plots for each cluster. Default is TRUE.
#' @param Heatmap Logical, if TRUE the function will generate a heatmap. Default is FALSE.
#' @param csv Logical, if TRUE the function will write CSV files for the top 10 markers and all markers per cluster. Default is FALSE.
#' @param pdf.height A numeric value specifying the height of the PDF. Default is 10.
#' @param pdf.width A numeric value specifying the width of the PDF. Default is 12.
#'
#' @return The function does not return anything. It saves the PDF and CSV files in the "output" directory.
#'
#' @seealso \code{\link[Seurat]{DimPlot}}, \code{\link[Seurat]{VlnPlot}}, \code{\link[Seurat]{DoHeatmap}}
#'
#' @importFrom grDevices pdf dev.off
#' @importFrom utils write.csv
#' @importFrom Seurat DimPlot PercentageFeatureSet VlnPlot FindAllMarkers DoHeatmap
#' @importFrom dplyr %>%, group_by, slice_max, top_n
#'
#' @export
#'


save_qc_file <- function(SeuFile, version = "v1.0.0", qc.plot = TRUE, Heatmap = FALSE, csv = FALSE, pdf.height = 10, pdf.width = 12){

  Filename <- as.character(substitute(SeuFile))

  if (!file.exists("output")) {
    dir.create("output")
  }

  pdf(paste0("output/",Filename, "_qc_file_", version, ".pdf", sep=""), width = pdf.width, height = pdf.height)

  # UMAP
  DefaultAssay(SeuFile) <- "RNA"
  umap <- DimPlot(SeuFile, reduction="umap", label=T, pt.size = 0.5)
  print(umap)

  if(qc.plot == TRUE){

    # VlnPlot : QC parameters of each cluster
    DefaultAssay(SeuFile) <- "RNA"
    SeuFile[["percent.mt"]] <- PercentageFeatureSet(SeuFile, pattern = "^MT-")
    qc <- VlnPlot(SeuFile, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    print(qc)

  }

  if(Heatmap == TRUE){

    # Heatmap
    markers <- FindAllMarkers(SeuFile, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

    markers %>%
      group_by(cluster) %>%
      slice_max(n = 10, order_by = avg_log2FC)

    markers %>%
      group_by(cluster) %>%
      top_n(n = 10, wt = avg_log2FC) -> top10

    hm <- DoHeatmap(SeuFile, features = top10$gene)
    print(hm)

  }

  dev.off()


  if(csv == TRUE){
    write.csv(top10, file = paste0("output/", Filename, "_Top10MarkersPerCluster_", version, ".csv"))
    write.csv(markers, file = paste0("output/", Filename, "_MarkersPerCluster_", version, ".csv"))
  }


  print(umap)

  if(qc.plot == TRUE){
    print(qc)
  }

  if(Heatmap == TRUE){
    print(hm)
  }

}
