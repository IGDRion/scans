#!/usr/bin/env Rscript

# Header ---------------------------------------------------------------------------------------------
## Aim : Format orthology files from orthoFinder
## Author(s): - Aurore BESSON
## Update : 

# Libraries -----------------------------------------------------------------------------------
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(tidyverse)))

# Variables -----------------------------------------------------------------------------------
## parse command line arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript format_orthology.R config.txt working_dir orthology_dir", call.=FALSE)
}
config_file <- args[1]
working_dir <- args[2]
orthology_dir <- args[3]


## set working directory
setwd(working_dir)

## create directories to store data
INPUTDIR <- "work/input_data"
dir.create(file.path(INPUTDIR, "orthology"), showWarnings = FALSE, recursive = TRUE)

# Input --------------------------------------------------------------------------------------

input_data <- read.table(config_file, header = TRUE, sep = ",", stringsAsFactors = FALSE)

ortho.list <- list.files(orthology_dir, ".tsv", full.names = T)

# Functions ---------------------------------------------------------------------------------------

# function to check if column containing `orthogroup` exists
contains_orthogroup <- function(x) grepl("orthogroup", x, ignore.case = TRUE)

# function to determine homology type
determine_relationship <- function(geneID_sp1, geneID_sp2) {
  sp1_count <- length(strsplit(geneID_sp1, ", ")[[1]])
  sp2_count <- length(strsplit(geneID_sp2, ", ")[[1]])
  
  if (sp1_count > 1 && sp2_count == 1) {
    return("ortholog_many2one")
  } else if (sp1_count == 1 && sp2_count > 1) {
    return("ortholog_one2many")
  } else if (sp1_count == 1 && sp2_count == 1) {
    return("ortholog_one2one")
  } else if (sp1_count > 1 && sp2_count > 1) {
    return("ortholog_many2many")
  } else {
    return("NA")
  }
}

# Script ---------------------------------------------------------------------------------------

## for each species, modify orthology files
for (ortho.path in ortho.list){
    # retrieve species names
    cat("orthology file in progress: ", basename(ortho.path), "\n")
    splitName <- strsplit(basename(ortho.path), "__v__|\\.")[[1]]
    query.name <- splitName[1]
    target.name <- splitName[2]
    
    # read file
    orthoDF <- as.data.frame(read.delim(ortho.path, header = T, stringsAsFactors = F))

    # check if `orthogroup` columns exist
    orthogroup_cols <- which(sapply(colnames(orthoDF), contains_orthogroup))
    ortho_col_names <- paste(colnames(orthoDF)[orthogroup_cols], collapse = ", ")
    if (length(orthogroup_cols) > 0) {
        # remove columns with `orthogroup`
        orthoDF <- orthoDF[, -orthogroup_cols]
        print(paste("Removal of columns containing 'orthogroup':", ortho_col_names))
    } else {
        print("No column containing 'orthogroup' in orthology file.")
    }

    # check if orthology file contains only 2 columns
    if (is.data.frame(orthoDF) && ncol(orthoDF) == 2) {
        # Rename the columns
        colnames(orthoDF) <- c("query", "target")
        print("Orthology columns renamed successfully.")
    } else {
        stop("Error: Orthology file does not have exactly 2 columns. Check your orthology file format.")
    }

    # create homology attribute
    orthoDF <- orthoDF %>%
        mutate(homolog_orthology_type = mapply(determine_relationship, 
            query, target)) %>%
        separate_rows(query, sep = ",") %>%
        separate_rows(target, sep = ",")

    # rename columns
    col2_name <- paste0(target.name,"_homolog")
    col3_name <- paste0(target.name,"_homolog_orthology_type")
    colnames(orthoDF) <- c("gene_id",col2_name,col3_name)

    # order dataframe according query gene_id
    orthoDF <- orthoDF %>%
        arrange(gene_id)

    # write modified orthology file
    file_name <- paste0(INPUTDIR,"/orthology/",query.name,"_",target.name,"_homology.tsv") 
    write.table(orthoDF, file_name, sep = "\t", row.names = FALSE, quote = FALSE)
    cat("Formatted orthology file created : ", file_name, "\n")
}   


    
