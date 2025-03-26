#!/usr/bin/env Rscript

# Header ---------------------------------------------------------------------------------------------
## Aim : Synteny analysis using method 2
## Author(s): - Fabien DEGALEZ
## Update : A. Besson

# Libraries -----------------------------------------------------------------------------------
suppressWarnings(suppressMessages(library(tidyr)))
suppressWarnings(suppressMessages(library(pbapply)))
suppressWarnings(suppressMessages(library(stringi)))
suppressWarnings(suppressMessages(library(stringr)))
suppressWarnings(suppressMessages(library(dplyr)))

# Variables -----------------------------------------------------------------------------------
## parse command line arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript perform_synteny_method2.R config.txt working_dir orthology_dir", call.=FALSE)
}
config_path <- args[1]
working_dir <- args[2]
orthology_dir <- args[3]

## set working directory
setwd(working_dir)

## create directories to store data in `method2``
dir_ortho <- paste0(working_dir,"/work/method2/orthoFeelnc")
dir_synt_byPair <- paste0(working_dir,"/work/method2/syntenyByPairFeelnc")
dir_synt_bySpecies <- paste0(working_dir,"/work/method2/mergedSyntenyBySpeciesFeelnc")

dir.create(dir_ortho, recursive = T)
dir.create(dir_synt_byPair, recursive = T)
dir.create(dir_synt_bySpecies, recursive = T)

# Input --------------------------------------------------------------------------------------

## configuration file
config <- read.table(config_path, header = TRUE, sep = ",", stringsAsFactors = FALSE)

# Functions ---------------------------------------------------------------------------------------

## function to extract target lncRNA linked to the homologous PCG
orthologous.fct <- function(i){
    PCG_lnc_1 <- lnc_1[lnc_1$feelLncPcgGnId == i, ]
    
    # Extraction of the lncRNA linked to the homologous PCG in species 2 
    homologous_i <- unlist(homology.file[(homology.file[, 1] == i | homology.file[, 2] == i), 1:2])
    homologous_i <- setdiff(homologous_i, i)
    if (!(identical(homologous_i, character(0)))) {
        PCG_lnc_2 <- lnc_2[lnc_2$feelLncPcgGnId == homologous_i, ]
    } 
    
    PCG_lnc_1_sub <- PCG_lnc_1[, c(1,3:5,7)]
    PCG_lnc_2_sub <- PCG_lnc_2[, c(1,3:5,7)]
    
    res_i <- merge(PCG_lnc_1_sub, PCG_lnc_2_sub,
                   by.x = "feelLncPcgClassName", by.y = "feelLncPcgClassName")
    
    if (nrow(res_i) == 0){
        return(NA)
    }
    
    res_i <- res_i[ , c(4,2,1,3,5,8,6,1,7,9)]
    return(res_i)
}

## function to determine class of orthology
whatOrthologies <- function(x){
    nbLncRNA_species1 <- str_count(x[1], ";")+1
    nbLncRNA_species2 <- str_count(x[2], ";")+1
    if (nbLncRNA_species1 == 1 & nbLncRNA_species2 == 1){
        return("one_to_one")
    } else if (nbLncRNA_species1 == 1 & nbLncRNA_species2 > 1) {
        return("one_to_many")
    } else if (nbLncRNA_species1 > 1 & nbLncRNA_species2 == 1) {
        return("many_to_one")
    } else if (nbLncRNA_species1 > 1 & nbLncRNA_species2 > 1) {
        if (nbLncRNA_species1 == nbLncRNA_species2) {
            return("many_to_many_strict")
        } else {
            return("many_to_many")
        }
    }
}

## function to remove duplicates
remove_duplicate <- function(x){
    toKeep <- !duplicated(as.character(str_split(x[1],";", simplify = T)))
    
    x[1] <- paste0(str_split(x[1],";", simplify = T)[toKeep], collapse = ";")
    x[2] <- paste0(str_split(x[2],";", simplify = T)[toKeep], collapse = ";")
    x[3] <- paste0(str_split(x[3],";", simplify = T)[toKeep], collapse = ";")
    return(c(x[1], x[2], x[3]))
}

## function to determine lnc configuration compared to his nearest PCG (position, strand, overlapping or not)
custom_config_1 <- function(x, thr_distance){
    #print(x)
    config_lnc <- x["feelLncPcgClassName"]
    
    strand_lnc <- x["strand"]
    end_lnc <- as.numeric(x["end"])
    start_lnc <- as.numeric(x["start"])
    
    
    strand_pcg <- x["feelLncPcg_strand"]
    end_pcg <- as.numeric(x["feelLncPcg_end"])
    start_pcg <- as.numeric(x["feelLncPcg_start"])
    
    
    if (config_lnc == "lincDivg"){
        distance <- min(abs(end_lnc - start_pcg), 
                        abs(start_lnc - end_pcg))
        if (distance <= thr_distance){
            return("lncgDivg")
        } else {
            return("lincDivg")
        }
    }
    
    if (config_lnc == "lincConv"){
        distance <- min(abs(end_lnc - start_pcg), 
                        abs(start_lnc - end_pcg))
        if (distance <= thr_distance){
            return("lncgConv")
        } else {
            return("lincConv")
        }
    }
    
    if (grepl("lincSS",config_lnc)){
        if (strand_lnc == "+"){
            if (end_lnc <= start_pcg){
                if ((start_pcg - end_lnc) <= thr_distance){
                    return("lncgSS.up")
                } else {
                    return("lincSS.up")
                }
            } else if (start_lnc >= end_pcg){
                if ((start_lnc - end_pcg ) <= thr_distance){
                    return("lncgSS.dw")
                } else {
                    return("lincSS.dw")
                }
            }
        } else if  (strand_lnc == "-"){
            if (start_lnc > end_pcg){
                if ((start_lnc - end_pcg) <= thr_distance){
                    return("lncgSS.up")
                } else {
                    return("lincSS.up")
                }
            } else if (start_pcg >= end_lnc){
                if ((start_pcg - end_lnc) <= thr_distance){
                    return("lncgSS.dw")
                } else {
                    return("lincSS.dw")
                }
            }
        } 
    }
    
    if (grepl("lncgSS",config_lnc)){
        if (strand_lnc == "+"){
            if (start_lnc < start_pcg){
                return("lncgSS.up")
            } else if (start_lnc >= start_pcg){
                return("lncgSS.dw")
            }
        } else if  (strand_lnc == "-"){
            if (start_lnc > start_pcg){
                return("lncgSS.up")
            } else if (start_lnc <= start_pcg){
                return("lncgSS.dw")
            }
            
            
        } 
    }
    
    if (grepl("lncgAS",config_lnc)){
        if (strand_lnc == "+"){
            if (start_lnc <= end_pcg & start_lnc >= start_pcg){
                return("lncgDivg")
            } else {
                return("lncgConv")
            }
        } else if  (strand_lnc == "-"){
            if (end_lnc <= end_pcg & end_lnc >= start_pcg){
                return("lncgDivg")
            } else {
                return("lncgConv")
            }
        } 
    }
    
    
    
    return("ERROR")
}

# function to uncollapse lnc
uncollapse_lnc <- function(x){
    x <- data.frame(t(x))
    x <- x[rep(1, str_count(x[1,1], ";")+1),]
    x[,1] <- as.character(str_split(x[1,1], ";", simplify = T))
    return(x)
}

# Script ---------------------------------------------------------------------------------------

## STEP1: Determine all possible comparisons of synteny
cat("STEP1: Determine synteny analysis to perform.", "\n")

toWork <- as.data.frame(expand_grid(config$completeName, config$completeName))
colnames(toWork) <- c("completeName_source", "completeName_target")
toWork$shortName_source <- config$shortName[match(toWork$completeName_source, config$completeName)]
toWork$shortName_target <- config$shortName[match(toWork$completeName_target, config$completeName)]
toWork <- toWork[toWork$shortName_source != toWork$shortName_target, ]

cat("Pair comparisons to perform by synteny method 2: ", "\n")
print(toWork)

## STEP2: identify orthologous lnc by pair
cat("STEP2: identification of orthologous lncRNA by pair.", "\n")

for (i in 1:nrow(toWork)){
    completeName_source <- as.character(toWork[i,1])
    completeName_target <- as.character(toWork[i,2])
    shortName_source <- as.character(toWork[i,3])
    shortName_target <- as.character(toWork[i,4])
    cat(shortName_source,
        "-", shortName_target, "\n")
    
    ## Import Homology file
    homology.path <- paste0(orthology_dir,"/", shortName_source,
                            "_", shortName_target, "_homology.tsv")
    cat("Homology file imported: ", homology.path, "\n")

    homology.file <- read.delim(homology.path, header = T, stringsAsFactors = F)
    homology.file <- homology.file[homology.file[,3] == "ortholog_one2one", ]
    
    ## Import Annotation file
    annot_1 <- read.delim(paste0("work/input_data/gnInfo/", completeName_source,
                                 "_gnInfo.tsv"),
                          header = T, stringsAsFactors = F)
    
    annot_2 <- read.delim(paste0("work/input_data/gnInfo/", completeName_target,
                                 "_gnInfo.tsv"),
                          header = T, stringsAsFactors = F)
    
    ## Import lncRNA file
    # lncRNA configuration files
    lnc_1 <- read.delim(paste0("work/method2/lncClassification/", completeName_source, "_lncConfiguration_feelncclassifier.tsv"),
                        header = T, stringsAsFactors = F)
    toKeep_lnc1 <- colnames(lnc_1)
    
    lnc_2 <- read.delim(paste0("work/method2/lncClassification/", completeName_target, "_lncConfiguration_feelncclassifier.tsv"),
                        header = T, stringsAsFactors = F)
    toKeep_lnc2 <- colnames(lnc_2)
    
    ## lncRNA selection
    lnc_1 <- lnc_1[!is.na(lnc_1$feelLncPcgClassName) 
                   & !grepl("unclassified", lnc_1$feelLncPcgClassName), ]
    
    lnc_2 <- lnc_2[!is.na(lnc_2$feelLncPcgClassName) 
                   & !grepl("unclassified", lnc_2$feelLncPcgClassName), ]
    
    ## Add annotation information
    lnc_1 <- merge(lnc_1, annot_1, by.x = "gnId", by.y = "gene_id", all.x = T )
    colnames(annot_1) <- paste0("feelLncPcg", "_", colnames(annot_1))
    lnc_1 <- merge(lnc_1, annot_1, by.x = "feelLncPcgGnId", by.y = "feelLncPcg_gene_id", all.x = T )
    
    
    lnc_2 <- merge(lnc_2, annot_2, by.x = "gnId", by.y = "gene_id", all.x = T )
    colnames(annot_2) <- paste0("feelLncPcg", "_", colnames(annot_2))
    lnc_2 <- merge(lnc_2, annot_2, by.x = "feelLncPcgGnId", by.y = "feelLncPcg_gene_id", all.x = T )
    
    
    # Configuration changement
    lnc_1$feelLncPcgClassName <- str_split(lnc_1$feelLncPcgClassName, "_", simplify = T)[,1]
    lnc_1$feelLncPcgClassName[lnc_1$feelLncPcgClassName == "lincSS"] <- "lncgSS"
    
    lnc_2$feelLncPcgClassName <- str_split(lnc_2$feelLncPcgClassName, "_", simplify = T)[,1]
    lnc_2$feelLncPcgClassName[lnc_2$feelLncPcgClassName == "lincSS"] <- "lncgSS"
    
    lnc_1$feelLncPcgClassName <- apply(lnc_1, 1, custom_config_1, 5000, simplify = T)
    lnc_2$feelLncPcgClassName <- apply(lnc_2, 1, custom_config_1, 5000, simplify = T)

    lnc_1 <- lnc_1[, toKeep_lnc1]
    lnc_2 <- lnc_2[, toKeep_lnc2]
    
    ### Creation of the final file 
    colNames <- c(paste0("PCG.", completeName_source),
                  paste0("lncRNA.", completeName_source),
                  paste0("conf.", completeName_source),
                  paste0("class.", completeName_source),
                  paste0("distance.", completeName_source),
                  paste0("PCG.", completeName_target),
                  paste0("lncRNA.", completeName_target),
                  paste0("conf.", completeName_target),
                  paste0("class.", completeName_target),
                  paste0("distance.", completeName_target))
    
    
    toExtract <- lnc_1$feelLncPcgGnId %in% homology.file[, 1] | lnc_1$feelLncPcgGnId %in% homology.file[, 2]
    PCG_homologous <- unique(lnc_1[toExtract, "feelLncPcgGnId"])
       
    res <- pbsapply(PCG_homologous[], orthologous.fct, USE.NAMES = T)
    res <- res[!is.na(res)]
    res <- do.call(rbind, res)
    
    colnames(res) <- colNames
    rownames(res) <- 1:nrow(res)
    filename <- paste0(dir_ortho,"/", shortName_source, "-", shortName_target, "_lncConfigurationHomology.tsv")
    write.table(res, filename, quote = F, row.names = F, sep = "\t", col.names = T)
    cat("Orthologous lncRNA saved in: ", filename, "\n")
    
    # aggregate results by pair
    res_agg <- aggregate(cbind(res[,2], res[,4], res[,5],
                               res[,7], res[,9], res[,10]),
                         list(res[,1], res[,3], res[,6], res[,8]),
                         paste, collapse = ";")
    
    colnames(res_agg) <- colnames(res[, c(1,3,6,8,2,4,5,7,9,10)])
    
    res_agg[,c(8,9,10)] <- data.frame(t(pbapply(res_agg[,c(8,9,10)], 1, remove_duplicate)), stringsAsFactors = F)
    res_agg[,c(5,6,7)] <- data.frame(t(pbapply(res_agg[,c(5,6,7)], 1, remove_duplicate)), stringsAsFactors = F)
    
    res_agg[,  paste0("orthology_type.", shortName_source, "-",shortName_target )] <- pbapply(res_agg[, c(5,8)], 1 , whatOrthologies)
    
    filename_agg <- paste0(dir_synt_byPair,"/", shortName_source, "-", shortName_target, "_lncConfigurationHomologyAggregated.tsv")
    write.table(res_agg, filename_agg, quote = F, row.names = F, sep = "\t", col.names = T)
    cat("Synteny by pair is saved in:", filename_agg, "\n")
}

# STEP3: Merge synteny analysis by species
cat("STEP3: Combine synteny analysis by query target.", "\n")

feelnc.list <- list.files(dir_synt_byPair, full.names = T, pattern = "_lncConfigurationHomologyAggregated.tsv")

for (i in 1:nrow(config)){
        
    shortName <- config$shortName[i]
    completeName <- config$completeName[i]
        
    cat(completeName,"\n")
        
    feelnc_interest.list <- feelnc.list[grep(paste0("/", shortName), feelnc.list)]
        
    res <- data.frame(matrix(NA, nrow = 0, ncol = 1), stringsAsFactors = F)
    colnames(res) <- paste0("lncRNA.", completeName)
        
    for (feelnc.path in feelnc_interest.list){
        feelnc.file <- read.delim(feelnc.path, header = T, stringsAsFactors = F, check.names = F)
        feelnc.file <- feelnc.file[, c(5,8,11)]
        test <- pbapply(feelnc.file, 1, uncollapse_lnc)
        test <- do.call(rbind.data.frame, test)
        res <- merge(res, test, 
                     by.x = paste0("lncRNA.", completeName),
                     by.y = paste0("lncRNA.", completeName), all = T)
    }
        
    tmp <- res[,grepl("orthology.type", colnames(res))]
    
    totalIndication <- function(x){
        return(sum(!is.na(x)))
    }
        
    profilIndication <- function(x){
        if (all(is.na(x))){
            return(NA)
        }
        x <- x[!is.na(x)]
        x <- data.frame(table(x), stringsAsFactors = F)
        x <- x[order(x$Freq),]            
        return(paste0(apply(x, 1, paste0, collapse=":"), collapse = ";"))
    }
    
    nbTotalIndication <- pbapply(tmp, 1, totalIndication)
    profilIndication <- pbapply(tmp, 1,  profilIndication)
        
    res$nbTotalIndication <- nbTotalIndication
    res$profilIndication <- profilIndication

    filename_synt <- paste0(dir_synt_bySpecies,"/",completeName,"_feelncMerged.tsv")
    write.table(res, filename_synt, quote = F, sep = '\t', row.names = F)
    cat("Merged synteny by species is saved in: ", filename_synt, "\n")
} 
