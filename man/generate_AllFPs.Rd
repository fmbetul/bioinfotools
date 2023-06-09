% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_allFPs.R
\name{generate_AllFPs}
\alias{generate_AllFPs}
\title{Generate Feature Plots for given genes across different clusters}
\usage{
generate_AllFPs(SeuFile, gene_list)
}
\arguments{
\item{SeuFile}{A Seurat object that contains the gene expression data.}

\item{gene_list}{A list of genes for which the feature plots need to be generated.}
}
\value{
Generates a PDF file containing the feature plots for the provided genes.
}
\description{
This function generates feature plots for each gene provided in the gene_list parameter across
different clusters in a Seurat object. It arranges these plots and saves them as a PDF file.
}
