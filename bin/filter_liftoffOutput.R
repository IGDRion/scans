#!/usr/bin/env Rscript

# Header ---------------------------------------------------------------------------------------------
## Aim : Filter liftoff output according cutoff (coverage) + keep only one copy for further analysis
## Author(s): - A. Besson
## Update : 

# Libraries -----------------------------------------------------------------------------------
suppressMessages(suppressWarnings(library(rtracklayer)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(stringr)))
suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(ggplot2)))

# Variables -----------------------------------------------------------------------------------
## parse command line arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 5) {
  stop("Usage: Rscript filter_liftoffOutput.R liftoff_out liftoff_modif coverage_cutoff identity_cutoff working_dir", call.=FALSE)
}
liftoff_out <- args[1]
liftoff_modif <- args[2]
coverage_cutoff <- args[3]
identity_cutoff <- args[4]
working_dir <- args[5]

## set working directory -> voir pour créer un dossier specifique pour tout enregistrer ?
setwd(working_dir)

# Input --------------------------------------------------------------------------------------

cat("Liftoff output file: ", liftoff_out, "\n")
cat("Coverage cut-off to apply: ", coverage_cutoff, "\n")
cat("Sequence identity cut-off to apply: ", identity_cutoff, "\n")

## apply cut-off on output from liftoff
liftoff_df <- as.data.frame(rtracklayer::import(liftoff_out))

## extract only genes
liftoff_genes <- liftoff_df %>% subset(type == "gene")

## visualization of coverage and sequence ID metrics
covPlot <- ggplot(liftoff_genes, aes(x = as.numeric(sequence_ID), y = as.numeric(coverage))) +
    geom_point() +
    theme_bw() +
    scale_x_continuous(name="Sequence identity") +
    scale_y_continuous(name="Coverage")

outplot <- paste0(working_dir,"/liftoff_plot_coverage_seqID.png")
cat("Plot of coverage in function of sequence identity saved in :",outplot, "\n")
ggsave(plot = covPlot, filename = outplot)

## remove extra-copies for following analysis if any
cat("Number of genes having extra copy (0 = no extra copy):", "\n")
table(liftoff_genes$extra_copy_number, useNA = "ifany")

## keep genes according cut-off (coverage, identity)
liftoff_genes_sub <- liftoff_genes %>%
    subset(as.numeric(coverage) >= coverage_cutoff & as.numeric(sequence_ID) >= identity_cutoff)
cat("Number of genes passing filters (coverage + sequence identity): ", length(liftoff_genes_sub$gene_id), "\n")

### filter df and save it as gtf
liftoff_cutoff <- liftoff_df %>%
    subset(gene_id %in% liftoff_genes_sub$gene_id) %>%
    subset(extra_copy_number == 0)

modified_ranges <- makeGRangesFromDataFrame(liftoff_cutoff, keep.extra.columns = TRUE)
rtracklayer::export(modified_ranges, liftoff_modif, format="gtf")

cat("Liftoff output file modified according cut-off: ", liftoff_modif, "\n")

