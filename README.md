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

### Generate Scaled Feature Plots

This function generates feature plots for a list of genes in multiple Seurat objects, scaling the plots by the maximum value of each gene across all objects. The generated plots are then saved as PDF files in an output directory.

```R
generate_scaled_FPs(seurat_objects, genes, filename, separate = TRUE)
```

- `seurat_objects`: A named list of Seurat objects.
- `genes`: A character vector of gene names for which to generate feature plots.
- `filename`: The filename (without extension) to use for the output PDF files.
- `separate`: A boolean indicating whether to save a separate PDF for each Seurat object. Default is `TRUE`.


### Generate Violin Plots for Given Genes Across Different Clusters

This function generates violin plots for each gene provided in the `gene_list` parameter across different clusters in a Seurat object. It arranges these plots and saves them as a PDF file.

```R
generate_allVlns(SeuFile, gene_list)
```

- `SeuFile`: A Seurat object that contains the gene expression data.
- `gene_list`: A list of genes for which the violin plots need to be generated.


### Plot Specific Cells from an Integrated Seurat Object

This function generates a UMAP plot for specific cells in a given sample that belong to a certain cluster from an integrated Seurat object.

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
