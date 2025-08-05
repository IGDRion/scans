#!/usr/bin/env Rscript

# Header ---------------------------------------------------------------------------------------------
## Aim : Generate global results output file for sequence alignment analysis (method 3)
## Author(s): - A. Besson
## Update : 

# Libraries -----------------------------------------------------------------------------------
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(tidyverse)))
suppressWarnings(suppressMessages(library(pbapply)))
suppressWarnings(suppressMessages(library(stringi)))
suppressWarnings(suppressMessages(library(stringr)))
suppressWarnings(suppressMessages(library(tidyr)))

# Variables -----------------------------------------------------------------------------------
## parse command line arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript getGlobalResultMethod3.R config.txt working_dir", call.=FALSE)
}
config <- args[1]
working_dir <- args[2]

## set working directory 
setwd(working_dir)

# Functions ---------------------------------------------------------------------------------------
## Function to uncollapse and count when multi-values separated by ";"
uncollapse_lnc <- function(x){
    x <- data.frame(t(x))
    x <- x[rep(1, str_count(x[1,1], ",")+1),]
    x[,1] <- as.character(str_split(x[1,1], ",", simplify = T))
    return(x)
}

## Function to count non-NA values in a row
totalIndication <- function(x){
    return(sum(!is.na(x)))
}

## Function to create a profile of orthology types
profilIndication <- function(x){
    if (all(is.na(x))){
        return(NA)
    }
    x <- x[!is.na(x)]
    x <- data.frame(table(x), stringsAsFactors = F)
    x <- x[order(x$Freq),]
    return(paste0(apply(x, 1, paste0, collapse=":"), collapse = ","))
}

# Input --------------------------------------------------------------------------------------
## all sequence alignment result files
seq.list <- list.files("scans_results/method3", "mapped_knownGenes.txt", full.names = T, recursive=TRUE)
seq.list.lnc <- seq.list[grep("lncRNA", seq.list)] # keep only lncRNA results

## config file
input_file <- read.table(config, header = TRUE, sep = ",", stringsAsFactors = FALSE)

# Generate global results file --------------------------------------------------------------------------------------
out_dir="scans_results/method3/mergedseqAligBySpecies"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

i <- 1
for (i in 1:nrow(input_file)){
    
    ## Import seqAlig data for each species 
    shortName <- input_file$shortName[i]
    completeName <- input_file$completeName[i]
    
    cat(completeName,"\n")
    
    seqAlig_interest.list <- grep(paste0("/", shortName, "_"), seq.list.lnc, value = TRUE)    
      
    ## Process each seqAlig file to obtain merged results by species
    res <- data.frame(matrix(NA, nrow = 0, ncol = 1), stringsAsFactors = F)
    colnames(res) <- paste0("lncRNA.", shortName)
    
    for (seqAlig.path in seqAlig_interest.list){
        speciesTarget <- strsplit(basename(seqAlig.path), "_")[[1]][3]
        seqAlig.file <- read.delim(seqAlig.path, header = T, stringsAsFactors = F)
        # format and uncollapse data
        seqAlig.file <- seqAlig.file[, c(1,6,8,15)]
        #seqAlig.file <- seqAlig.file %>% subset(target_gene_biotype %in% lncRNA.regex)
        colnames(seqAlig.file) <- c(paste0("lncRNA.", shortName), paste0("geneID.",speciesTarget), paste0("biotype.",speciesTarget), paste0("orthology_type.",shortName,"-",speciesTarget))
        test <- pbapply(seqAlig.file, 1, uncollapse_lnc)
        test <- do.call(rbind.data.frame, test)
        # merge formatted seqAlig data by species
        res <- merge(res, test, 
                     by.x = paste0("lncRNA.", shortName),
                     by.y = paste0("lncRNA.", shortName), all = T)
    }
    
    tmp <- res[,grepl("orthology.type", colnames(res))]
    tmp2 <- res[,grepl("biotype", colnames(res))]
    
    
    # add nb of orthologies + orthology class for each gene by species
    nbTotalIndication <- pbapply(tmp, 1, totalIndication)
    profilInd <- pbapply(tmp, 1,  profilIndication)
    profilBiotype <- pbapply(tmp2, 1,  profilIndication)
    res$nbTotalIndication <- nbTotalIndication
    res$profilIndication <- profilInd
    res$biotypeIndication <- profilBiotype
    
    # save merged results
    filename <- paste0(out_dir,"/",completeName,"_seqAligMerged.tsv")
    write.table(res, filename, quote = F, sep = '\t', row.names = F)
    cat("Generation of global results analysis for method 3: ", filename, "\n")
}