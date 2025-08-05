#!/usr/bin/env Rscript

# Header ---------------------------------------------------------------------------------------------
## Aim : Generate global results output file for one species in the 3 methods
## Author(s): - A. Besson
## Update : 

# Libraries -----------------------------------------------------------------------------------
suppressMessages(suppressWarnings(library(tidyr)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(stringr)))
suppressWarnings(suppressMessages(library(pbapply)))

# Variables -----------------------------------------------------------------------------------
## parse command line arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript getGlobalSummary.R shortNameSpecies working_dir", call.=FALSE)
}
species <- args[1]
working_dir <- args[2]

## set working directory 
setwd(working_dir)

# Functions ---------------------------------------------------------------------------------------
# function to get max value by converting values in numeric or NA if non numerical values
# by taking in account mulitple values separated by ";"
safe_max <- function(x) {
  nums <- suppressWarnings(as.numeric(trimws(x)))
  if (all(is.na(nums))) return(NA)
  max(nums, na.rm = TRUE)
}

## Function to uncollapse and count when multi-values separated by ";"
uncollapse_lnc <- function(x){
    x <- data.frame(t(x))
    x <- x[rep(1, str_count(x[1,1], ";")+1),]
    x[,1] <- as.character(str_split(x[1,1], ";", simplify = T))
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

# Data ---------------------------------------------------------------------------------------
## directory containing all method detailed results
meth1_dir <- "method1/syntenyByPair"
meth2_dir <- "method2/syntenyByPairFeelnc"
meth3_dir <- paste0("method3/",species)

## all result files
meth1.list <- list.files(meth1_dir, paste0("^",species), full.names = T, recursive=TRUE)

meth2.list <- list.files(meth2_dir, paste0(species,"-"), full.names = T, recursive=TRUE)

meth3.list <- list.files(meth3_dir, "mapped_knownGenes.txt", full.names = T, recursive=TRUE)
meth3.list <- meth3.list[grep("lncRNA", meth3.list)] # keep only lncRNA results

# Global summary by pair ---------------------------------------------------------------------------------------
## directory containing all global results by pair
summaryByPair <- paste0("results_summary/summaryByPair/",species)
dir.create(summaryByPair, showWarnings = FALSE, recursive = TRUE)

## write summary results for each pair
for (meth1.file in meth1.list){
    # pair analyzed
    pairName <- str_replace(basename(meth1.file),"_synteny.tsv","")
    cat("\n","Pair analysis in progress: ", pairName, "\n")
    speciesQuery <- strsplit(basename(meth1.file), "_")[[1]][1]
    speciesTarget <- strsplit(basename(meth1.file), "_")[[1]][2]
    
    # METHOD 1
    
    cat("\n","Result file Method 1: ", meth1.file, "\n")
    ## read data
    meth1_data <- read.delim(meth1.file, header = T, stringsAsFactors = F)
    meth1_data <- meth1_data[, c(1,8,9,10,17,18,20,21)]
    ## separate lncRNA query (one per line)
    meth1_data <- meth1_data %>%
        separate_rows(c(1,2,3), sep = ";") %>%
        mutate(across(c(2,3), as.numeric))
    ## get max value for distance lncRNA-PCG
    meth1_data <- meth1_data %>%
        mutate(
            max_col5 = sapply(strsplit(as.character(.[[5]]), ";"), safe_max),
            max_col6 = sapply(strsplit(as.character(.[[6]]), ";"), safe_max),
            lncRNA_PCG_max_distance_query = pmax(.[[2]],.[[3]], na.rm = TRUE),
            lncRNA_PCG_max_distance_target = pmax(max_col5, max_col6, na.rm = TRUE)) %>%
        select(-max_col5, -max_col6, -c(2,3,5,6))
    ## reorder & rename
    meth1_data <- meth1_data %>% select(c(1,5,2,6,4,3))
    colnames(meth1_data) <- c(paste0("lncRNA.", speciesQuery), 
                              paste0("M1_lnc_PCG_max_distance.",speciesQuery),
                              paste0("M1_lncRNA.",speciesTarget),
                              paste0("M1_lnc_PCG_max_distance.",speciesTarget),
                              paste0("M1_orientation.",speciesQuery,".",speciesTarget),
                              paste0("M1_orthology.",speciesQuery,".",speciesTarget))
    ## remove orthology with zero
    meth1_data <- meth1_data %>% 
        filter(!str_detect(.[[6]], 'zero'))
    
    # METHOD 2
    
    ## extract species name
    pairName2 <- str_replace(basename(pairName),"_","-")
    meth2.file <- list.files(meth2_dir, pattern=pairName2, full.names=TRUE)
    cat("\n","Result file Method 2: ", meth2.file, "\n")
    ## read data
    meth2_data <- read.delim(meth2.file, header = T, stringsAsFactors = F)
    meth2_data <- meth2_data[, c(5,2,7,8,10,11)]
    ## separate lncRNA query (one per line)
    meth2_data <- meth2_data %>%
        separate_rows(c(1,3), sep = ";") %>%
        mutate(across(c(3), as.numeric))
    ## get max value for distance lncRNA-PCG
    meth2_data <- meth2_data %>%
        mutate(max_col5 = sapply(strsplit(as.character(.[[5]]), ";"), safe_max)) %>%
        select(c(1,2,3,4,7,6))

    colnames(meth2_data) <- c(paste0("lncRNA.", speciesQuery), 
                          paste0("M2_FEELnc_class.", speciesQuery),
                          paste0("M2_lnc_PCG_distance.",speciesQuery),
                          paste0("M2_lncRNA.",speciesTarget),
                          paste0("M2_lnc_PCG_distance.",speciesTarget),
                          paste0("M2_orthology.",speciesQuery,".",speciesTarget))
    
    # METHOD 3
    
    pairName3 <- str_replace(basename(pairName),"_","_to_")
    meth3.file <- meth3.list[grep(pairName3, meth3.list)]
    cat("\n","Result file Method 3: ", meth3.file, "\n")
    ## read data
    meth3_data <- read.delim(meth3.file, header = T, stringsAsFactors = F)
    meth3_data <- meth3_data[, c(1,5,6,8,10,11,15)]
    colnames(meth3_data) <- c(paste0("lncRNA.", speciesQuery), 
                          paste0("M3_frac_overlap_max.", speciesQuery),
                          paste0("M3_geneID.",speciesTarget),
                          paste0("M3_biotype.",speciesTarget),
                          paste0("M3_frac_overlap_max.",speciesTarget),
                          paste0("M3_strand_match.",speciesQuery,".",speciesTarget),
                          paste0("M3_orthology.",speciesQuery,".",speciesTarget))
    
    # MERGE RESULTS + ADD CODE
    
    ## merge results from the 3 methods
    M1_M2 <- merge(meth1_data, meth2_data, by.x = names(meth1_data)[1],
                   by.y = names(meth2_data)[1], all = T)
    all_df <- merge(M1_M2, meth3_data, by.x = names(M1_M2)[1],
                   by.y = names(meth3_data)[1],all = T)
    ## add lnc class code
    all_df <- all_df %>% 
    mutate(lnc_class = case_when(
        # code 1: M1|M2|M3
        (!(is.na(.[[6]])) & is.na(.[[11]]) & is.na(.[[17]])) |
        (!(is.na(.[[11]])) & is.na(.[[6]]) & is.na(.[[17]])) | 
        (!(is.na(.[[17]])) & is.na(.[[11]]) & is.na(.[[6]])) ~ 1,
        # Code 2 : M1 + M2
        (!(is.na(.[[6]])) & !(is.na(.[[11]])) & is.na(.[[17]])) ~ 2,
        # Code 3 : M2 + M3
        (!(is.na(.[[17]])) & !(is.na(.[[11]])) & is.na(.[[6]])) ~ 3,
        # Code 4 : M1 + M3
        (!(is.na(.[[17]])) & !(is.na(.[[6]])) & is.na(.[[11]])) ~ 4,
        # Code 5 : M1 + M2 + M3
        (!(is.na(.[[17]])) & !(is.na(.[[6]])) & !(is.na(.[[11]]))) ~ 5,
        # no orthology
        TRUE ~ 0))
    ## add lnc corresponding id to check if the lnc id is similar between methods
    all_df <- all_df %>%
        mutate(lnc_corresp = case_when(
        # Code 1 : M1 | M2 | M3
        lnc_class == 1 ~ "NA",
        # Code 2 : M1 + M2
        (lnc_class == 2 & lengths(mapply(intersect, strsplit(.[[3]], ","), 
                                        strsplit(.[[9]], ","))) > 0) ~ "same_lncRNA",
        # Code 3 : M2 + M3
        (lnc_class == 3 & lengths(mapply(intersect, strsplit(.[[9]], ","), 
                                         strsplit(.[[13]], ","))) > 0) ~ "same_lncRNA",
        # Code 4 : M1 + M3
        (lnc_class == 4 & lengths(mapply(intersect, strsplit(.[[3]], ","), 
                                        strsplit(.[[13]], ","))) > 0) ~ "same_lncRNA",
        # Code 5 : M1 + M2 + M3
        (lnc_class == 5 & lengths(mapply(intersect, strsplit(.[[3]], ","), strsplit(.[[9]], ","))) > 0 &
         lengths(mapply(intersect, strsplit(.[[3]], ","), strsplit(.[[13]], ","))) > 0) ~ "same_lncRNA",
    
        TRUE ~ "different_lncRNA"))
    ## rename columns
    names(all_df)[names(all_df) == 'lnc_class'] <- paste0("lnc_class.",speciesQuery,".",speciesTarget)
    names(all_df)[names(all_df) == 'lnc_corresp'] <- paste0("lnc_corresp.",speciesQuery,".",speciesTarget)
    ## write file
    fileName <- paste0(summaryByPair,"/",speciesQuery,"_",speciesTarget,"_globalResults.tsv")
    write.table(all_df, file = fileName, quote = F, row.names = F,  sep = '\t')
    cat("\n","Summary results in: ", fileName, "\n")
    print(table(all_df$lnc_class, useNA = "ifany"))
}


# Global summary for one species ---------------------------------------------------------------------------------------
all.list <- list.files(summaryByPair, "_globalResults.tsv", full.names = T, recursive=TRUE)

## create dataframe
res <- data.frame(matrix(NA, nrow = 0, ncol = 1), stringsAsFactors = F)
colnames(res) <- paste0("lncRNA.",species)
 
## for each by-pair summary
for (pair.path in all.list){
    cat("\n","File to analyzed:",basename(pair.path),"\n")
    speciesTarget <- strsplit(basename(pair.path), "_")[[1]][2]
    
    #read file
    pair.file <- read.delim(pair.path, header = T, stringsAsFactors = F)
        
    # select, format & uncollapse
    pair.file <- pair.file %>% select(c(1,6,11,14,17,18,19))
    pair.file[[4]] <- ifelse((pair.file[[4]] == "lncRNA,mRNA" | pair.file[[4]] == "mRNA,lncRNA"),"both",pair.file[[4]])
    pair.file[[6]] <- sub("^","class ",pair.file[[6]])
    pair.file$M1_species <- ifelse(is.na(pair.file[[2]]),NA,speciesTarget)
    pair.file$M2_species <- ifelse(is.na(pair.file[[3]]),NA,speciesTarget)
    pair.file$M3_species <- ifelse(is.na(pair.file[[5]]),NA,speciesTarget)
    names(pair.file)[names(pair.file) == 'M1_species'] <- paste0("M1_species.",speciesQuery,".",speciesTarget)
    names(pair.file)[names(pair.file) == 'M2_species'] <- paste0("M2_species.",speciesQuery,".",speciesTarget)
    names(pair.file)[names(pair.file) == 'M3_species'] <- paste0("M3_species.",speciesQuery,".",speciesTarget)
        
    test <- pbapply(pair.file, 1, uncollapse_lnc)
    test <- do.call(rbind.data.frame, test)
    
    # merge by query lncRNA id to have all target species infos in one file
    res <- merge(res, test, by.x = names(res)[1], by.y = names(test)[1], all = T)
}
## retrieve results by method
tmpM1sp <- res[,grepl("M1_species.", colnames(res))]
tmpM2sp <- res[,grepl("M2_species.", colnames(res))]
tmpM3sp <- res[,grepl("M3_species.", colnames(res))]
tmpM1 <- res[,grepl("M1_orthology.", colnames(res))]
tmpM2 <- res[,grepl("M2_orthology.", colnames(res))]
tmpM3 <- res[,grepl("M3_orthology.", colnames(res))]
tmpBioT <- res[,grepl("M3_biotype.", colnames(res))]
tmpClass <- res[,grepl("lnc_class.", colnames(res))]
tmpCorr <- res[,grepl("lnc_corresp.", colnames(res))]
    
    
## add nb of orthologies + lnc class for each gene by species
res$M1_nbTotalSpecies <- pbapply(tmpM1, 1, totalIndication)
res$M1_orthology <- pbapply(tmpM1, 1, profilIndication)

res$M2_nbTotalSpecies <- pbapply(tmpM2, 1, totalIndication)
res$M2_orthology <- pbapply(tmpM2, 1, profilIndication)
    
res$M3_nbTotalSpecies <- pbapply(tmpM3, 1, totalIndication)
res$M3_orthology <- pbapply(tmpM3, 1, profilIndication)
res$M3_biotype <- pbapply(tmpBioT, 1, profilIndication)

res$allMethods_lncClass <- pbapply(tmpClass, 1, profilIndication)
res$allMethods_lncCorresp <- pbapply(tmpCorr, 1, profilIndication)

## add species name if orthology for each gene by method
res <- res %>%
    mutate(M1_species = apply(select(., starts_with("M1_species")), 1, function(row) {
            vals <- row[!is.na(row) & row != "NA"]
            if(length(vals) == 0) NA else paste(vals, collapse = ";")}),
        M2_species = apply(select(., starts_with("M2_species")), 1, function(row) {
            vals <- row[!is.na(row) & row != "NA"]
            if(length(vals) == 0) NA else paste(vals, collapse = ";")}),
        M3_species = apply(select(., starts_with("M3_species")), 1, function(row) {
            vals <- row[!is.na(row) & row != "NA"]
            if(length(vals) == 0) NA else paste(vals, collapse = ";")}))
    
## select columns
final_df <- res %>% 
            select(1,M1_nbTotalSpecies,M1_orthology,M1_species,
                    M2_nbTotalSpecies,M2_orthology,M2_species,
                    M3_nbTotalSpecies,M3_orthology,M3_species,M3_biotype,
                    allMethods_lncClass,allMethods_lncCorresp)

## save merged results
filename <- paste0("results_summary/",species,"_conservation.tsv")
write.table(final_df, filename, quote = F, sep = '\t', row.names = F)
cat(species,"conservation results is saved in: ", filename, "\n")