#!/usr/bin/env Rscript

# Header ---------------------------------------------------------------------------------------------
## Aim : Format files to perform synteny analysis
## Author(s): - Fabien DEGALEZ
## Update : A. Besson

# Libraries -----------------------------------------------------------------------------------
suppressWarnings(suppressMessages(library(rtracklayer)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(tidyverse)))

# Variables -----------------------------------------------------------------------------------
## parse command line arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript format_gtf.R config.txt working_dir", call.=FALSE)
}
config_file <- args[1]
working_dir <- args[2]

## set working directory
setwd(working_dir)

## create directories to store data
INPUTDIR <- "work/input_data"
dir.create(file.path(INPUTDIR, "gnInfo"), showWarnings = FALSE, recursive = TRUE)

# Input --------------------------------------------------------------------------------------

input_file <- read.table(config_file, header = TRUE, sep = ",", stringsAsFactors = FALSE)
input_file

# Script ---------------------------------------------------------------------------------------

## read each gtf and keep only genes
gnInfoList = list()

for (i in 1:nrow(input_file)) {
    GTF.name <- input_file$completeName[i]
    cat("GTF modification in progress: ", GTF.name, "\n")

    gtf_data <- as.data.frame(rtracklayer::import(input_file$pathToGTF[i]))

    # check if gene_biotype is present
    if (!"gene_biotype" %in% colnames(gtf_data)) {
        if ("gene_type" %in% colnames(gtf_data)) {
            # if gene_type, rename the column
            gtf_data <- gtf_data %>% rename(gene_biotype = gene_type)
            print("'gene_type' column has been renamed as 'gene_biotype'.")
        } else {
            print("Warning : No 'gene_biotype' nor 'gene_type' have been found. Check your GTF file.")
        }
    } else {
        print("GTF file format OK.")
    }

    # select genes and wanted features
    gtf_genes <- gtf_data %>%
        subset(type == "gene") %>%
        dplyr::select_if(names(.) %in% c("gene_id","gene_name","gene_biotype","seqnames","start","end","strand"))
    
    # add species name
    gtf_genes$species <- GTF.name

    # write output
    output_file <- paste0(INPUTDIR,"/gnInfo/",GTF.name, "_gnInfo.tsv")
    write.table(gtf_genes, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
  
    cat("Created file: ", output_file, "\n")

    # merge all species
    gnInfoList[[GTF.name]] <- gtf_genes
}

# bind all gtf files and write the output
gnInfo <- do.call(rbind, gnInfoList)
write.table(gnInfo, paste0(INPUTDIR,"/allMerged_gnInfo.tsv"), 
            quote = F, sep = '\t', row.names = F)

