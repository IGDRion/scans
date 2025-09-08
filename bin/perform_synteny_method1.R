#!/usr/bin/env Rscript

# Header ---------------------------------------------------------------------------------------------
## Aim : Synteny analysis using method 1
## Author(s): - Fabien DEGALEZ
## Update : A. Besson

# Libraries -----------------------------------------------------------------------------------
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(pbapply)))
suppressWarnings(suppressMessages(library(stringi)))
suppressWarnings(suppressMessages(library(stringr)))
suppressWarnings(suppressMessages(library(tidyr)))

# Variables -----------------------------------------------------------------------------------
## parse command line arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
  stop("Usage: perform_synteny_method1.R config.txt working_dir orthology_dir", call.=FALSE)
}
config_file <- args[1]
working_dir <- args[2]
orthology_dir <- args[3]

## set working directory
setwd(working_dir)

## create directories to store data in `method1``
SYNTDIR <- "scans_results/method1"
dir.create(file.path(SYNTDIR, "lncBetwPcg"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(SYNTDIR, "syntenyByPair"), showWarnings = FALSE)
dir.create(file.path(SYNTDIR, "mergedSyntenyBySpecies"), showWarnings = FALSE)

# Input --------------------------------------------------------------------------------------

## configuration file
input_file <- read.table(config_file, header = TRUE, sep = ",", stringsAsFactors = FALSE)

## gnInfo files
gnInfo.list <- list.files(file.path(working_dir,"scans_results/input_data/gnInfo"), "gnInfo.tsv", full.names = T)

## complete biotype labelling
PCG.regex <- c("protein_coding","mRNA")
lncRNA.regex <- c("lncRNA","lincRNA", "antisense", "sense_overlapping", "sense_intronic", "lnc_RNA")


# Functions ---------------------------------------------------------------------------------------

## Function to determine PCG id surrounding lnc
pcgIdEachSide <- function(y){
    x <- y["gene_id"]
    ## lncRNA position
    chrLNC <- y["seqnames"]
    startLNC <- as.numeric(y["start"])
    endLNC <- as.numeric(y["end"])
    strandLNC <- y["strand"]
    
    startLNC_real <- startLNC
    endLNC_real <- endLNC
    
    ## Left PCG
    leftPCG <- GTF[GTF$seqnames == chrLNC & GTF$start < startLNC & GTF$simpleBiotype == "pcg" , ]
    if (nrow(leftPCG) > 0) {
        leftPCG <- leftPCG[, c("seqnames", "start", "end", "strand", "gene_id")]
        leftPCG$diff <- startLNC_real - leftPCG$start
        leftPCG <- leftPCG[leftPCG$diff == min(leftPCG$diff), ]
        leftPCG <- leftPCG[1, ]
    } else {
        leftPCG <- data.frame(t(rep(NA, 6)), stringsAsFactors = F)
        colnames(leftPCG) <- c("seqnames", "start", "end", "strand", "gene_id", "diff")
    }
    
    ## Right PCG
    rightPCG <- GTF[GTF$seqnames == chrLNC & GTF$start >= startLNC & GTF$simpleBiotype == "pcg", ]
    if (nrow(rightPCG) > 0) {
        rightPCG <- rightPCG[, c("seqnames", "start", "end", "strand", "gene_id")]
        rightPCG$diff <- abs(startLNC_real - rightPCG$start)
        rightPCG <- rightPCG[rightPCG$diff == min(rightPCG$diff), ]
        rightPCG <- rightPCG[1, ]
    }else {
        rightPCG <- data.frame(t(rep(NA, 6)), stringsAsFactors = F)
        colnames(rightPCG) <- c("seqnames", "start", "end", "strand", "gene_id", "diff")
    }
    
    ## Configuration type of the lncRNA compared to the tw PCGs
    if (any(!is.na(leftPCG)) & any(!is.na(rightPCG))) {
        LNC_start <- startLNC
        LNC_end <- endLNC
        leftPCG_start <- leftPCG$start
        leftPCG_end <- leftPCG$end
        rightPCG_start <- rightPCG$start
        rightPCG_end <- rightPCG$end
        
        # Case 1 
        if  (((leftPCG_end >= LNC_start) & (LNC_start >= leftPCG_start)) &
             ((leftPCG_end >= LNC_end) & (LNC_end >= leftPCG_start))){
            LNC_conf <- "contained_PCG_Left"
            # Case 2
        } else if (((leftPCG_end >= LNC_start) & (LNC_start >= leftPCG_start)) &
                   ((leftPCG_end <= LNC_end) & (LNC_end <= rightPCG_start))){
            LNC_conf <- "overlap_PCG_Left"
            # Case 3
        } else if (((leftPCG_end >= LNC_start) & (LNC_start >= leftPCG_start)) &
                   ((rightPCG_end >= LNC_end) & (LNC_end >= rightPCG_start))){
            LNC_conf <- "overlap_PCG_both"
            # Case 4 
        } else if (((rightPCG_start >= LNC_start) & (LNC_start >= leftPCG_end)) &
                   ((rightPCG_start >= LNC_end) & (LNC_end >= leftPCG_end))){
            LNC_conf <- "between"
            # Case 5
        } else if (((rightPCG_start >= LNC_start) & (LNC_start >= leftPCG_end)) &
                   ((rightPCG_start <= LNC_end) & (LNC_end <= rightPCG_end))){
            LNC_conf <- "overlap_PCG_Right"
            # Case 6 
        } else if (((rightPCG_start <= LNC_start) & (LNC_start <= rightPCG_end)) &
                   ((rightPCG_start <= LNC_end) & (LNC_end <= rightPCG_end))){
            LNC_conf <- "contained_PCG_Right"
            # Case 7 
        } else if (((leftPCG_start <= LNC_start) & (LNC_start <= leftPCG_end)) &
                   ((LNC_end >= rightPCG_end))){
            LNC_conf <- "overlap_PCG_left_AND_PCG_right_contained"
            # Case 8
        } else if (((leftPCG_end <= LNC_start) & (LNC_start <= rightPCG_start)) &
                   ((LNC_end >= rightPCG_end))){
            LNC_conf <- "PCG_right_contained"
        } else {
            LNC_conf <- "Error/other"
        }
    } else {
        LNC_conf <- NA
    }

    return(c(x, leftPCG$gene_id, rightPCG$gene_id,
             strandLNC, leftPCG$strand, rightPCG$strand,
             LNC_conf, leftPCG$diff, rightPCG$diff))
}

## Function to check strands of orthologous PCG : same, reverse or discordant
test_PCG <- function(PCG_source, PCG_target){
    if (PCG_source == PCG_target){
        return("same")
    }
    PCG_source <- str_replace_all(PCG_source,"\\+","A")
    PCG_source <- str_replace_all(PCG_source,"\\-","\\+")
    PCG_source <- str_replace_all(PCG_source,"A","\\-")
    if (PCG_source == PCG_target){
        return("reverse")
    }
    return("discordant")
}

## Function to find orthologous lncRNAs in target species based on surrounding PCG
find_lncRNA_orth_byPCGcouple <- function(x){
    
    # Source PCG pair
    PCG_left <- x[1]
    PCG_right <- x[2]
    table_source_tmp <- tableSource[tableSource[,2] %in% PCG_left & tableSource[,3] %in% PCG_right, ] 
    colnames_table_source <- colnames(table_source_tmp)
    
    # Find orthologous target PCG pair
    PCG_orthologous_left <- homology[match(PCG_left, homology[,1]), 2]
    PCG_orthologous_right<- homology[match(PCG_right, homology[,1]), 2]
    table_target_tmp <- tableTarget[((tableTarget[, 2] %in% PCG_orthologous_left) | (tableTarget[, 3] %in% PCG_orthologous_left)) & 
                                        ((tableTarget[, 2] %in% PCG_orthologous_right) | (tableTarget[, 3] %in% PCG_orthologous_right)), ]
    colnames_table_target <- colnames(table_target_tmp)
    
    # Test concordance PCGs strand
    PCG_source <- paste0(c(unique(table_source_tmp[, 5]),unique(table_source_tmp[, 6])), collapse = "")
    PCG_target <- paste0(c(unique(table_target_tmp[, 5]),unique(table_target_tmp[, 6])), collapse = "")
    
    # Process source table data
    table_source_toPaste <- table_source_tmp[, c(1,4,7,8,9)]    
    table_source_toPaste_forward <- table_source_toPaste[table_source_toPaste[,2] == "+",]
    table_source_toPaste_reverse <- table_source_toPaste[table_source_toPaste[,2] == "-",] 
    
    # Split by strand to avoid Discordant
    if (nrow(table_source_toPaste_forward) > 0) {
        
        table_source_toPaste_forward <- t(apply(table_source_toPaste_forward, 2, function(x){return(paste0(x, collapse = ";"))}))
        table_source_toPaste_forward <- cbind(table_source_tmp[1, c(2,3,5,6)], table_source_toPaste_forward)
        table_source_toPaste_forward <- table_source_toPaste_forward[, c(5,1,2,6,3,4,7,8,9)]
        colnames(table_source_toPaste_forward) <- colnames_table_source
    } else {
        table_source_toPaste_forward <- data.frame(t(rep(NA, ncol(table_source_tmp))), stringsAsFactors = F)
        colnames(table_source_toPaste_forward) <- colnames_table_source
    }
    
    if (nrow(table_source_toPaste_reverse) > 0) {
        table_source_toPaste_reverse <- t(apply(table_source_toPaste_reverse, 2, function(x){return(paste0(x, collapse = ";"))}))
        table_source_toPaste_reverse <- cbind(table_source_tmp[1, c(2,3,5,6)], table_source_toPaste_reverse)
        table_source_toPaste_reverse <- table_source_toPaste_reverse[, c(5,1,2,6,3,4,7,8,9)]
        colnames(table_source_toPaste_reverse) <- colnames_table_source
    } else {
        table_source_toPaste_reverse <- data.frame(t(rep(NA, ncol(table_source_tmp))), stringsAsFactors = F)
        colnames(table_source_toPaste_reverse) <- colnames_table_source
    }
    ## Orthologous in the target species
    
    table_target_toPaste <- table_target_tmp[, c(1,4,7,8,9)]
    table_target_toPaste_forward <- table_target_toPaste[table_target_tmp[,4] == "+",]
    table_target_toPaste_reverse <- table_target_toPaste[table_target_tmp[,4] == "-",]
    
    if (nrow(table_target_toPaste_forward) > 0) {
        ## Concatenation
        table_target_toPaste_forward <- t(apply(table_target_toPaste_forward, 2, function(x){return(paste0(x, collapse = ";"))}))
        table_target_toPaste_forward <- cbind(table_target_tmp[1, c(2,3,5,6)], table_target_toPaste_forward)
        table_target_toPaste_forward <- table_target_toPaste_forward[, c(5,1,2,6,3,4,7,8,9)]
        colnames(table_target_toPaste_forward) <- colnames_table_target
    } else {
        table_target_toPaste_forward <- data.frame(t(rep(NA, ncol(table_target_tmp))), stringsAsFactors = F)
        colnames(table_target_toPaste_forward) <- colnames_table_target
    }
    
    if (nrow(table_target_toPaste_reverse) > 0) {
        table_target_toPaste_reverse <- t(apply(table_target_toPaste_reverse, 2, function(x){return(paste0(x, collapse = ";"))}))
        table_target_toPaste_reverse <- cbind(table_target_tmp[1, c(2,3,5,6)], table_target_toPaste_reverse)
        table_target_toPaste_reverse <- table_target_toPaste_reverse[, c(5,1,2,6,3,4,7,8,9)]
        colnames(table_target_toPaste_reverse) <- colnames_table_target
        
    } else {
        table_target_toPaste_reverse <- data.frame(t(rep(NA, ncol(table_target_tmp))), stringsAsFactors = F)
        colnames(table_target_toPaste_reverse) <- colnames_table_target
    }
    
    if (nrow(table_target_tmp) > 0){

        if (PCG_orthologous_left == unique(table_target_tmp[,3])){
            info <- "symmetry"
            PCG_target <- stri_reverse(PCG_target)
        } else {
            info <- "normal"
        }    
    }
    strand_couple <- test_PCG(PCG_source, PCG_target)
    
    if (strand_couple == "reverse"){
        res_forward <- cbind(table_source_toPaste_forward, table_target_toPaste_reverse)
        res_reverse <- cbind(table_source_toPaste_reverse, table_target_toPaste_forward)
    } else {
        res_forward <- cbind(table_source_toPaste_forward, table_target_toPaste_forward)
        res_reverse <- cbind(table_source_toPaste_reverse, table_target_toPaste_reverse)
    }
    
    res <- rbind(res_forward, res_reverse)
    res$PCG_strandCouple <- strand_couple
    return(res)
}

## Function to determine orthology categories: one/many_to_one/many/zero
orthology_categories <- function(x){
    nb_lncRNA_source <- str_count(x[1], ";") + 1
    nb_lncRNA_target <- str_count(x[10], ";") + 1
    
    if ((nb_lncRNA_source == 1) & (is.na(nb_lncRNA_target))){
        return("one_to_zero")
    } else if ((nb_lncRNA_source > 1) & (is.na(nb_lncRNA_target))){
        return("many_to_zero")
    } else if ((nb_lncRNA_source == 1) & (nb_lncRNA_target == 1)){
        return("one_to_one")
    } else if ((nb_lncRNA_source == 1) & (nb_lncRNA_target > 1)){
        return("one_to_many")
    } else if ((nb_lncRNA_source > 1) & (nb_lncRNA_target == 1)){
        return("many_to_one")
    } else if ((nb_lncRNA_source > 1) & (nb_lncRNA_target > 1)){
        if ((nb_lncRNA_source == nb_lncRNA_target)){
            return("many_to_many_sameNumber")
        } else {
            return("many_to_many_DifferentNumber")
        }
    } else {
        return("Error")
    }
}

## Function to collapse strand for the "many" cases
strand_collapsing <- function(x){
    strand <- as.character(str_split(x, ";", simplify = T))
    strand_res <- c(strand[1])
    initialStrand <- strand[1]
    for (i in strand) {
        if (i != initialStrand)
            strand_res <- c(strand_res, i)
        initialStrand <- i 
    }
    return(paste0(strand_res, collapse = ";"))
}

## Function to determine orientation class for orthologous pcg-lnc-pcg group
strand_orientation <- function(x){
    ## Extracxtion of orientation from source species 
    strand_lncRNA_source <- strand_collapsing(x[4])
    strand_source <- paste0(c(x[5], strand_lncRNA_source, x[6] ), collapse = "")
    ## Extracxtion of orientation from target species 
    strand_lncRNA_target <- strand_collapsing(x[13])
    strand_target <- paste0(c(x[14], strand_lncRNA_target, x[15] ), collapse = "")
    ## Look for a matching group
    resSource <- strand_case_possibility[match(strand_source, strand_case_possibility$combinaison),]
    resTarget <- strand_case_possibility[match(strand_target, strand_case_possibility$combinaison),]
    
    if (nchar(strand_source) != nchar(strand_target)){
        orientation_group <- NA
        orientation_class <- "discordant_multi"
    } else {
        if ((nchar(strand_source) > 3) &  (nchar(strand_target) > 3)){
            if (strand_source == strand_target){
                orientation_group <- NA
                orientation_class <- "same_multi"
            } else {
                strand_source <- str_replace_all(strand_source, "\\+", "a")
                strand_source <- str_replace_all(strand_source, "\\-", "\\+")
                strand_source <- str_replace_all(strand_source, "a", "\\-")
                if (strand_source == strand_target){
                    orientation_group <- NA
                    orientation_class <- "reverse_multi"
                } else {
                    orientation_group <- NA
                    orientation_class <- "discordant_multi"
                }
            }
        } else if  ((nchar(strand_source) == 3) &  (nchar(strand_target) == 3)){
            if (resSource$combinaison == resTarget$combinaison){
                orientation_group <- resSource$group
                orientation_class <- "same"
            } else if ((resSource$combinaison != resTarget$combinaison) & (resSource$group == resTarget$group)){
                orientation_group <- resSource$group
                orientation_class <- "reverse"
            } else {
                orientation_group <- NA
                orientation_class <- "discordant"
            }
        }
    }
    return(c(orientation_class, orientation_group))
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
    return(paste0(apply(x, 1, paste0, collapse=":"), collapse = ";"))
}

# Script ---------------------------------------------------------------------------------------

## STEP1: Determine all possible comparisons of synteny
cat("STEP1: Determine synteny analysis to perform.", "\n")

comparisons <- expand.grid(shortName_source = input_file$shortName,
                           shortName_target = input_file$shortName,
                           stringsAsFactors = FALSE) %>%
    # add completeName
    left_join(input_file, by = c("shortName_source" = "shortName")) %>%
    rename(completeName_source = completeName) %>%
    left_join(input_file, by = c("shortName_target" = "shortName")) %>%
    rename(completeName_target = completeName) %>%
    # filter 
    filter(shortName_source != shortName_target) %>%
    select(shortName_source, shortName_target, completeName_source, completeName_target)

cat("Pair comparisons to perform by synteny method 1: ", "\n")
print(comparisons)

## STEP2: find surrounding PCG for each lncRNA by species
cat("STEP2: find surrounding PCG for each lncRNA by species.", "\n")

for (gnInfo.path in gnInfo.list){
    # use gnInfo file
    GTF.name <- gsub("_gnInfo.tsv","",basename(gnInfo.path), fixed=T)
    cat(GTF.name, "\n")
    GTF <- read.delim(gnInfo.path, header = T, stringsAsFactors = F)

    # Selection of PCG and lncRNA only and format df
    GTF <- GTF %>% 
        subset(gene_biotype %in% c(PCG.regex,lncRNA.regex)) %>%
        #rename(seqname = seqnames) %>%
        dplyr::select_if(names(.) %in% c("seqnames","start","end","strand","gene_id","gene_name","gene_biotype"))
    cat("Total protein coding and lncRNA genes analyzed: ", length(GTF$gene_id), "\n")
    
    # add simpleBiotype
    GTF <- GTF %>% 
        mutate(simpleBiotype = case_when(
            gene_biotype %in% lncRNA.regex ~ "lnc",
            gene_biotype %in% PCG.regex ~ "pcg",
            TRUE ~ "na"))
    
    # select lnc only
    listLNC <- GTF %>% subset(gene_biotype %in% lncRNA.regex)
    cat("Total lncRNA genes analyzed: ", length(listLNC$gene_id), "\n")
    
    # table with lnc between PCG
    boundedLNC <- data.frame(t(pbapply(listLNC, 1, pcgIdEachSide)), stringsAsFactors = F)
    colnames(boundedLNC) <- c(paste0("lncRNA.", GTF.name), paste0("PCG_left.", GTF.name), paste0("PCG_right.", GTF.name),
                              paste0("lncRNA_strand.", GTF.name), paste0("PCG_left_strand.", GTF.name), paste0("PCG_right_strand.", GTF.name),
                              paste0("lncRNA_conf.", GTF.name), paste0("lncRNA_PCG_left_distance.", GTF.name), paste0("lncRNA_PCG_right_distance.", GTF.name))
    
    filename <- paste0(SYNTDIR, "/lncBetwPcg/",GTF.name,"_lncRNAbetweenPcg.tsv")
    write.table(boundedLNC, filename, quote = F, sep = '\t', row.names = F)
    cat("Table PCG-lnc-PCG saved in: ", filename, "\n")
}

# Table with possible strand configurations PCG-lnc-PCG
strand_case_possibility <- data.frame(c("+++",
                                        "++-", "+-+", "-++",
                                        "+--", "-+-", "--+",
                                        "---"),
                                      c(1,
                                        2, 3, 4,
                                        4, 3, 2,
                                        1), 
                                      stringsAsFactors = F)
colnames(strand_case_possibility) <- c("combinaison", "group")

# STEP3: perform synteny for each pair of species (2 anchors synteny: PCG-lnc_PCG)
cat("STEP3: Perform synteny for each species pair: 2 anchors based method (PCG-lnc-PCG)", "\n")

i <- 2
for (i in 1:nrow(comparisons)){
    cat(i,"/",nrow(comparisons), ": ")
    
    ## extract query/target species names for comparison
    shortName_source <- as.character(comparisons[i,1])
    shortName_target <- as.character(comparisons[i,2])
    completeName_source <- as.character(comparisons[i,3])
    completeName_target <- as.character(comparisons[i,4])
    
    cat(shortName_source,"-", shortName_target, "\n")
    
    ## Import data
    source.path <- paste0(SYNTDIR, "/lncBetwPcg/",completeName_source,"_lncRNAbetweenPcg.tsv")
    tableSource <- read.delim(source.path, header = T, stringsAsFactors = F)
    cat("query input data: ", source.path, "\n")
    
    target.path <- paste0(SYNTDIR, "/lncBetwPcg/", completeName_target,"_lncRNAbetweenPcg.tsv")
    tableTarget <- read.delim(target.path, header = T, stringsAsFactors = F)
    cat("target input data: ", target.path, "\n")
    
    homology.path <- paste0(orthology_dir,"/", shortName_source,
                            "_", shortName_target, "_homology.tsv")
    homology <- read.delim(homology.path, header = T, stringsAsFactors = F)
    cat("homology input data: ", homology.path, "\n")
    
    ## Filter and process homology data 
    # keep only one-to-one orthologs
    homology <- homology[homology[, 3] == "ortholog_one2one", ]
    
    # create correspondance table for source species
    correspondance <- data.frame(tableSource[, 1], tableSource[, 2] %in% homology[, 1] ,
                                 tableSource[, 3] %in% homology[, 1], stringsAsFactors = F)
    
    colnames(correspondance) <- c("lncRNA_source_id", "PCG_left_source_orthologous", "PCG_right_source_orthologous")
    
    # filter lncRNAs with orthologous flanking PCGs
    lncRNA_concerned <- correspondance[correspondance$PCG_left_source_orthologous == T &
                                           correspondance$PCG_right_source_orthologous == T   , "lncRNA_source_id"]
    tableSource <- tableSource[tableSource[, 1] %in% lncRNA_concerned, ]
    tableSource <- tableSource[tableSource[,2] != tableSource[,3],]
    
    tableSource_PCGcouple <- unique(tableSource[, c(2,3)])
    
    # find orthologous lncRNAs based on PCG couples
    orthologous_lncRNA_cases <- pbapply(tableSource_PCGcouple, 1, find_lncRNA_orth_byPCGcouple, simplify = T)
    
    # Process the results of orthologous lncRNA search
    orthologous_lncRNA_cases <- do.call(rbind, lapply(orthologous_lncRNA_cases, data.frame, row.names = NULL))
    orthologous_lncRNA_cases <- orthologous_lncRNA_cases[!apply(orthologous_lncRNA_cases, 1, function(x){all(is.na(x))}), ]
    orthologous_lncRNA_cases <- orthologous_lncRNA_cases[!is.na(orthologous_lncRNA_cases[,1]),]
    
    # Determine orthology type for each case
    orthologous_lncRNA_cases[,  paste0("orthology_type.", completeName_source, "-",completeName_target )] <- pbapply(orthologous_lncRNA_cases, 1, orthology_categories)
    
    # Split data based on orthology type (no orthology vs at least one orthology)
    tableSource_orientation_strand <- orthologous_lncRNA_cases[(orthologous_lncRNA_cases$orthology_type != "one_to_zero") & (orthologous_lncRNA_cases$orthology_type != "many_to_zero"), ]
    tableSource_no_orientation_strand <- orthologous_lncRNA_cases[(orthologous_lncRNA_cases$orthology_type == "one_to_zero") | (orthologous_lncRNA_cases$orthology_type == "many_to_zero"), ]
    
    # Determine the class and group orientation
    toAdd <- data.frame(t(pbapply(tableSource_orientation_strand, 1, strand_orientation)), stringsAsFactors = F)
    colnames(toAdd) <- c(paste0("orientation_class.", completeName_source, "-",completeName_target ),
                         paste0("orientation_group.", completeName_source, "-",completeName_target ))
    tableSource_orientation_strand <- cbind(tableSource_orientation_strand, toAdd)
    
    # Add NA values for orientation in cases without orthologous lncRNAs
    tableSource_no_orientation_strand[, paste0("orientation_class.", completeName_source, "-",completeName_target )] <- NA
    tableSource_no_orientation_strand[ ,paste0("orientation_group.", completeName_source, "-",completeName_target )] <- NA
    
    ## Merge the processed data and save
    tableSource <- rbind(tableSource_orientation_strand, tableSource_no_orientation_strand)
    
    filename <- paste0(SYNTDIR, "/syntenyByPair/",shortName_source,"_",shortName_target,"_synteny.tsv")
    write.table(tableSource, filename, quote = F, sep = '\t', row.names = F)
    cat("Synteny by pair is saved in: ", filename, "\n")
}

# STEP4: Merge synteny analysis by species
cat("STEP4: Combine synteny analysis by query target.", "\n")

synteny.list <- list.files(paste0(SYNTDIR, "/syntenyByPair"), full.names = T, pattern = "synteny.tsv")

i <- 1
for (i in 1:nrow(input_file)){
    
    ## Import synteny data for each species 
    shortName <- input_file$shortName[i]
    completeName <- input_file$completeName[i]
    
    cat(completeName,"\n")
    
    synteny_interest.list <- grep(paste0("/", shortName, "_"), synteny.list, value = TRUE)    
    print(synteny_interest.list)
    
    ## Process each synteny file to obtain merged results by species
    res <- data.frame(matrix(NA, nrow = 0, ncol = 1), stringsAsFactors = F)
    colnames(res) <- paste0("lncRNA.", completeName)
    
    for (synteny.path in synteny_interest.list){
        synteny.file <- read.delim(synteny.path, header = T, stringsAsFactors = F)
        # format and uncollapse data
        synteny.file <- synteny.file[, c(1,10,20)]
        test <- pbapply(synteny.file, 1, uncollapse_lnc)
        test <- do.call(rbind.data.frame, test)
        # merge formatted synteny data by species
        res <- merge(res, test, 
                     by.x = paste0("lncRNA.", completeName),
                     by.y = paste0("lncRNA.", completeName), all = T)
    }
    
    # format for non ortholog genes
    res[res == "one_to_zero"] <- NA
    res[res == "many_to_zero"] <- NA
    
    tmp <- res[,grepl("orthology.type", colnames(res))]
    
    
    # add nb of orthologies + orthology class for each gene by species
    nbTotalIndication <- pbapply(tmp, 1, totalIndication)
    profilInd <- pbapply(tmp, 1,  profilIndication)
    res$nbTotalIndication <- nbTotalIndication
    res$profilIndication <- profilInd
    
    # save merged results
    filename <- paste0(SYNTDIR, "/mergedSyntenyBySpecies/",completeName,"_syntenyMerged.tsv")
    write.table(res, filename, quote = F, sep = '\t', row.names = F)
    cat("Merged synteny by species is saved in: ", filename, "\n")
}