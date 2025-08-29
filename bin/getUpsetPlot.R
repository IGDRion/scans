#!/usr/bin/env Rscript

# Header ---------------------------------------------------------------------------------------------
## Aim : Generate upset plot by species pair (methods comparison) or by method (one species vs all others)
## Author(s): - A. Besson
## Update : 

# Libraries -----------------------------------------------------------------------------------
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(stringr)))
suppressMessages(suppressWarnings(library(UpSetR)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(ggpubr)))

# Variables -----------------------------------------------------------------------------------
## parse command line arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 5) {
  stop("Usage: Rscript getUpsetPlot.R upsetMode shortNameSpecies specificArg working_dir config.txt", call.=FALSE)
}
upsetMode <- args[1]
species <- args[2]
specificArg <- args[3]
working_dir <- args[4]
config <- args[5]

## set working directory + create directory to write figures
setwd(working_dir)
dir.create("plots", showWarnings = FALSE, recursive = TRUE)

# Functions ---------------------------------------------------------------------------------------

# function to generate upset plot by pair
# function to generate upset plot by pair
getUpsetPlotByPair <- function(mergedFilePath, method, speciesCN){
    # read file
    merged_df <- as.data.frame(read.delim(mergedFilePath, header = T, stringsAsFactors = F))
    # according method retrieve list of one-to-one orthologous query lncRNA id
    if (method %in% c("met1","met2")){
        sub_df <- merged_df %>% filter(!str_detect(profilIndication, "many"))

        # columns 'lncRNA' to select
        lncRNA_cols <- grep("^lncRNA\\.", names(sub_df), value = TRUE)
        lncRNA_cols <- setdiff(lncRNA_cols, paste0("lncRNA.",speciesCN))
        
        # extract species name + ID lnc ref if orthologs
        listOrthoGenes <- purrr::map(lncRNA_cols, function(col) {
            sub_df[[1]][!is.na(sub_df[[col]])]})
        
        # name according suffixe behind "lncRNA."
        names(listOrthoGenes) <- sub("lncRNA\\.", "", lncRNA_cols)

    } else if (method == "met3"){
        sub_df <- merged_df %>% 
            filter(!str_detect(profilIndication, "many")) %>%
            filter(str_detect(biotypeIndication, "lncRNA"))
        lncRNA_cols <- grep("^geneID\\.", names(sub_df), value = TRUE)
        
        # extract species name + ID lnc ref if orthologs
        listOrthoGenes <- purrr::map(lncRNA_cols, function(col) {
            sub_df[[1]][!is.na(sub_df[[col]])]})
        
        # name according suffixe behind "lncRNA."
        names(listOrthoGenes) <- sub("geneID\\.", "", lncRNA_cols)

     }else {
        cat("Wrong method. Use met1|met2|met3")
    }

    # upset plot
    return(listOrthoGenes)
}

# Data ---------------------------------------------------------------------------------------
# config file to retrieve species name
config_df <- read.table(config, sep = ',', h = T)

# retrieve complete name
speciesCN <- config_df$completeName[config_df$shortName == species]
speciesCN

# Upset plot  ---------------------------------------------------------------------------------------
if (upsetMode == "byMethod"){
    ## data
    mergedSynt1.list <- list.files("method1/mergedSyntenyBySpecies", speciesCN, full.names = T)
    mergedSynt2.list <- list.files("method2/mergedSyntenyBySpeciesFeelnc", speciesCN, full.names = T)
    mergedSeqA.list <- list.files("method3/mergedseqAligBySpecies", speciesCN, full.names = T)
    cat("\n","Input files: ", "\n")
    cat(mergedSynt1.list,"\n")
    cat(mergedSynt2.list,"\n")
    cat(mergedSeqA.list,"\n")

    ## generate upset plot
    upsetList <- getUpsetPlotByPair(mergedSynt1.list, specificArg, speciesCN)

    filename <- paste0("plots/upsetPlot_",species,"_ortho1to1lnc_",specificArg,".png")

    png(filename, width = 2000, height = 1600, res = 200)
    print(upset(fromList(upsetList), keep.order = TRUE, 
      text.scale = 2, point.size = 4, line.size = 2, nsets = 8))
    dev.off()
    cat("\n","Upset plot by method saved: ", filename, "\n")

} else if (upsetMode == "byPair"){
    ## data
    pair <- paste0(species,"_",specificArg)
    pair2 <- paste0(species,"-",specificArg)
    pair3 <- paste0(species,"_to_",specificArg,"_mapped_knownGenes.txt")

    synt.list <- list.files("method1/syntenyByPair", pair, full.names = T)
    synt2.list <- list.files("method2/syntenyByPairFeelnc", pair2, full.names = T)
    seq.list <- list.files("method3", pair3, full.names = T, recursive=TRUE)
    seq.list <- seq.list[grep("lncRNA", seq.list)]

    cat("\n","Input files: ", "\n")
    cat(synt.list,"\n")
    cat(synt2.list,"\n")
    cat(seq.list,"\n")

    ## generate upset plot
    ### synteny df M1
    syntByPair <- as.data.frame(read.delim(synt.list, header = T, stringsAsFactors = F))
    syntByPair <- syntByPair %>% select(c(1,8,9,10,17,18,19,20))
    colnames(syntByPair) <- c("lnc_query","query_distance_PCG_left_M1","query_distance_PCG_right_M1",
                            "lnc_target_M1","target_distance_PCG_left_M1","target_distance_PCG_right_M1",
                            "class_strand_M1","class_orthology_M1")

    ### synteny df M2
    syntByPair2 <- as.data.frame(read.delim(synt2.list, header = T, stringsAsFactors = F))
    syntByPair2 <- syntByPair2 %>% select(c(5,2,7,8,10,11))
    colnames(syntByPair2) <- c("lnc_query","lnc_query_conf","query_distance_PCG_M2",
                            "lnc_target_M2","target_distance_PCG_M2","class_orthology_M2")

    # seq align df
    seqAlignByPair <- as.data.frame(read.delim(seq.list, header = T, stringsAsFactors = F))
    seqAlignByPair <- seqAlignByPair %>% select(c(1,5,6,8,10,11,15))
    colnames(seqAlignByPair) <- c("lnc_query","query_frac_overlap_M3","gene_target_M3","gene_target_biotype_M3","target_frac_overlap_M3","class_strand_M3","class_orthology_M3")

    ### subset to keep only one-to-one orthology with lnc target
    seqAlignByPair_sub <- seqAlignByPair %>%
        subset(class_orthology_M3 == "one_to_one" & gene_target_biotype_M3 == "lncRNA")

    syntByPair_sub <- syntByPair %>% 
        subset(class_orthology_M1 == "one_to_one")

    syntByPair2_sub <- syntByPair2 %>% 
        subset(class_orthology_M2 == "one_to_one")

    ### merge synteny and seq alignment results
    mergedDF <- merge(x = syntByPair_sub, y = seqAlignByPair_sub,
                    by = "lnc_query", all = TRUE)

    mergedDF2 <- merge(x = mergedDF, y = syntByPair2_sub,
                    by = "lnc_query", all = TRUE)
    ### columns 'class_orthology' to select
    method_cols <- grep("^class_orthology_", names(mergedDF2), value = TRUE)

    ### extract method + ID lnc ref if orthologs
    listOrthoGenesByM <- purrr::map(method_cols, function(col) {
        # retrieve lnc id if orthology = one_to_one (not NA)
        mergedDF2$lnc_query[!is.na(mergedDF2[[col]])]})

    ### rename column names
    names(listOrthoGenesByM) <- sub("class_orthology_", "", method_cols)
    renaming <- c("M1" = "method 1", "M2" = "method 2", "M3" = "method 3")
    names(listOrthoGenesByM) <- renaming[names(listOrthoGenesByM)]

    ### save plot
    filename <- paste0("plots/upsetPlot_",species,"_vs_",specificArg,"_ortho1to1lnc_allMethods.png")

    png(filename, width = 2000, height = 1600, res = 200)
    print(upset(fromList(listOrthoGenesByM), keep.order = T, 
                sets = c("method 3", "method 2", "method 1"),
                  text.scale=2, point.size = 4, line.size = 2))
    dev.off()

    cat("\n","Upset plot by pair saved: ", filename, "\n")
}

