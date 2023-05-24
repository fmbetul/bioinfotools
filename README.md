# bioinfotools

The `bioinfotools` package provides essential bioinformatics tools for the analysis and visualization of single-cell RNA sequencing data.

## Installation

You can install the `bioinfotools` package from GitHub using the `devtools` package:

``` r
devtools::install_github("fmbetul/bioinfotools")
```

## Usage

To use the package, load it into your R session:

``` r
library(bioinfotools)
```

For more information on available functions and their usage, please refer to the package documentation.


### 1. Generate Scaled Feature Plots

This function generates feature plots for a list of genes in multiple Seurat objects, scaling the plots by the maximum value of each gene across all objects. The generated plots are then saved as PDF files in an output directory.

```R
generate_scaled_FPs(seurat_objects, genes, filename, separate = TRUE)
```

- `seurat_objects`: A named list of Seurat objects.
- `genes`: A character vector of gene names for which to generate feature plots.
- `filename`: The filename (without extension) to use for the output PDF files.
- `separate`: A boolean indicating whether to save a separate PDF for each Seurat object. Default is `TRUE`.


### 2. Generate Feature Plots for Given Genes

This function generates feature plots for each gene provided in the `gene_list` parameter in a Seurat object. It arranges these plots and saves them as a PDF file.

```R
generate_AllFPs(SeuFile, gene_list)
```

- `SeuFile`: A Seurat object that contains the gene expression data.
- `gene_list`: A list of genes for which the feature plots need to be generated.


### 3. Generate Violin Plots for Given Genes Across Different Clusters

This function generates violin plots for each gene provided in the `gene_list` parameter across different clusters in a Seurat object. It arranges these plots and saves them as a PDF file.

```R
generate_allVlns(SeuFile, gene_list)
```

- `SeuFile`: A Seurat object that contains the gene expression data.
- `gene_list`: A list of genes for which the violin plots need to be generated.


### 4. Save QC File

This function saves UMAP plots and quality control (QC) files, with options to include QC plots and heatmaps, to the output folder.

```R
save_qc_file(SeuFile, version = "v1.0.0", qc.plot = TRUE, Heatmap = FALSE, csv = FALSE, pdf.height = 10, pdf.width = 12)
```

- `SeuFile`: A Seurat object that contains the data.
- `version`: A string that provides the version number for the QC file. Default is "v1.0.0".
- `qc.plot`: A logical indicating whether to generate QC plots for each cluster. Default is `TRUE`.
- `Heatmap`: A logical indicating whether to generate a heatmap. Default is `FALSE`.
- `csv`: A logical indicating whether to write CSV files for the top 10

 markers and all markers per cluster. Default is `FALSE`.
- `pdf.height`: A numeric value specifying the height of the PDF. Default is 10.
- `pdf.width`: A numeric value specifying the width of the PDF. Default is 12.


### 5. Plot Cells of a Specified Sample and Cluster Type in Integrated Seurat Object

This function plots cells of a specified sample and cluster type using a UMAP plot.

```R
plot_cells_from_sample(sample, ClusterName, integrated)
```

- `sample`: A Seurat object that contains the sample data.
- `ClusterName`: A string that specifies the cell type or cluster name to be plotted.
- `integrated`: An integrated Seurat object that contains the merged data from multiple samples.


### 6. Plot Specific Cells from an Integrated Seurat Object

The `plot_cells_from_integrated` function generates a UMAP plot for specific cells in a given sample that belong to a certain cluster from an integrated Seurat object.

```R
plot_cells_from_integrated(integrated, ClusterName, sample)
```

- `integrated`: A Seurat object. This is an integrated or merged Seurat object that includes data from multiple samples.
- `ClusterName`: A character string specifying the name of the cell type or cluster whose cells are to be visualized.
- `sample`: A Seurat object that contains scRNA-seq data for a specific sample.



## Dependencies

The `bioinfotools` package depends on the following packages:

- Seurat
- dplyr
- ggplot2
- ggpubr
- magrittr
- grDevices
- utils
- rlang

## License

The `bioinfotools` package is licensed under the GPL-3 license.

## Issues and Contributions

If you encounter any issues with the `bioinfotools` package or would like to contribute to its development, please [open an issue](https://github.com/fmbetul/bioinfotools/issues) on GitHub.

## Contact

For any further inquiries or questions, please feel free to contact the package maintainer, F. M. Betul Erol, at fmab_erol@hotmail.com.
