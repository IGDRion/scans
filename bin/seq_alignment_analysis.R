#!/usr/bin/env Rscript

# Header ---------------------------------------------------------------------------------------------
## Aim : Sequence alignment analysis (method 3)
## Author(s): - A. Besson
## Update : 

# Libraries -----------------------------------------------------------------------------------
suppressMessages(suppressWarnings(library(rtracklayer)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(stringr)))

# Variables -----------------------------------------------------------------------------------
## parse command line arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 5) {
  stop("Usage: Rscript seq_alignment_analysis.R target.gtf query.gtf overlap_query_to_target.bed unmapped_features.txt working_dir", call.=FALSE)
}
query_gtf <- args[1]
target_gtf <- args[2]
overlap_bed <- args[3]
unmap_txt <- args[4]
working_dir <- args[5]

## set working directory -> voir pour créer un dossier specifique pour tout enregistrer ?
setwd(working_dir)

# Functions ---------------------------------------------------------------------------------------

# function to extract species info from gtf
getGTFinfo <- function(gtf_path, species){
    gtf_df <- as.data.frame(rtracklayer::import(gtf_path))
    # check if transcript_biotype is present
    if (!"transcript_biotype" %in% colnames(gtf_df)) {
        if ("transcript_type" %in% colnames(gtf_df)) {
            # if `type`, rename the column by `biotype`
            gtf_df <- gtf_df %>% 
                rename(gene_biotype = gene_type) %>% 
                rename(transcript_biotype = transcript_type)
            print("'gene/transcript_type' column has been renamed as 'transcript/gene_biotype'.")
        } else {
            print("Warning : No 'transcript/gene_biotype' nor 'transcript/gene_type' have been found. Check your GTF file.")
        }
    } else {
        print("GTF file format OK.")
    }
    # check if gene_name is present or create column with NAs
    if (!"gene_name" %in% colnames(gtf_df)) {
        gtf_df <- gtf_df %>% add_column(gene_name = NA)
        print("'gene_name' column has been created.")
    } else {
        print("'gene_name' already in GTF.")
    }
    # get tx length by summing exon lengths
    exons_df <- gtf_df %>%
        subset(type == "exon") %>%
        group_by(transcript_id) %>%
        summarise(transcript_length = sum(width))
    # merge to have all info
    info_df <- gtf_df %>%
        subset(type == "transcript") %>%
        dplyr::select(transcript_id, transcript_biotype, gene_name, gene_id, gene_biotype) %>%
        left_join(exons_df, by = "transcript_id")
    colnames(info_df) <- c(paste0(species,"_transcript_id"),
                           paste0(species,"_transcript_biotype"),
                           paste0(species,"_gene_name"),
                           paste0(species,"_gene_id"),
                           paste0(species,"_gene_biotype"),
                           paste0(species,"_transcript_length"))
    
    return(info_df)
}

# Input --------------------------------------------------------------------------------------

species1 <- str_extract(basename(overlap_bed), "(?<=overlap_).*?(?=_to_)")
species2 <- str_extract(basename(overlap_bed), "(?<=_to_).*?(?=\\.bed)")

## retrieve gtf info for both species 
cat("Get transcript/gene info from: ", query_gtf, "\n")
query_df <- getGTFinfo(query_gtf, "query")

cat("Get transcript/gene info from: ", target_gtf, "\n")
target_df <- getGTFinfo(target_gtf, "target")

## format intersection output file from bedtools
cat("Format bedtools intersection file: ", overlap_bed, "\n")
overlap_df <- read.table(overlap_bed)
overlap_df <- overlap_df %>% 
    dplyr::select(V1, V2, V3, V4, V6, V10,
                  V13, V14, V15, V16, V18, V22, 
                    V25)
colnames(overlap_df) <- c("seqnames_liftoff","start_liftoff","end_liftoff","query_transcript_id","strand_liftoff","nb_exons_liftoff",
                            "seqnames_annot", "start_annot","end_annot","target_transcript_id","strand_annot","nb_exons_annot",
                            "overlap_length")


# not aligned genes - liftoff output ----------------------------------------------------------

cat("STEP1 - Query genes unmapped on target genes : ","\n")

cat("liftoff output file containing unmapped query genes : ", unmap_txt,"\n")
unmap_df <- read.table(unmap_txt, h = F)
colnames(unmap_df) <- "query_gene_id"

## add gtf info (to change if species name is used instead of query/target !!!)
query_genes <- query_df %>%
    select(query_gene_name, query_gene_id, query_gene_biotype) %>%
    distinct()
unmap_df <- unmap_df %>%
    left_join(query_genes, by = "query_gene_id")
print(table(unmap_df$query_gene_biotype, useNA = "ifany"))

out_file_zero <- paste0(species1,"_to_",species2,"_unmapped_genes.txt")
write.table(x = unmap_df, file = out_file_zero, quote = F, row.names = F, sep= "\t")
cat("File containing query genes unmapped on target genes : ", out_file_zero,"\n")

# transcript analysis -------------------------------------------------------------------------

cat("STEP 2 - Intersection analysis at transcript level", "\n")
## transcripts aligned but not overlapping known genes in annotation
wo_overlap <- overlap_df %>% 
    subset(seqnames_annot == ".") %>% 
    select(seqnames_liftoff, start_liftoff, end_liftoff, query_transcript_id, 
           strand_liftoff, nb_exons_liftoff) %>%
    left_join(query_df, by = "query_transcript_id")

out_file_wo <- paste0(species1,"_to_",species2,"_mapped_unknownTranscripts.txt")
write.table(x = wo_overlap, file = out_file_wo, quote = F, row.names = F, sep= "\t")
cat("File with",species1, "transcripts mapped on",species2, "but not overlapping with known transcripts in annotation: ", out_file_wo, "\n")

## transcripts aligned and overlapping known genes in annotation
with_overlap <- overlap_df %>% 
    subset(seqnames_annot != ".") %>%
    left_join(query_df, by = "query_transcript_id" ) %>%
    left_join(target_df, by = "target_transcript_id")

## add info if overlap in the same/different strand
with_overlap$strand_match <- ifelse(
    with_overlap$strand_liftoff == with_overlap$strand_annot,
    "same_strand", "antisense"
)

out_file_with <- paste0(species1,"_to_",species2,"_mapped_knownTranscripts.txt")
write.table(x = with_overlap, file = out_file_with, quote = F, row.names = F, sep= "\t")
cat("File with",species1, "transcripts mapped on",species2, "and overlapping with known transcripts in annotation: ", out_file_with, "\n")
   
# gene analysis -------------------------------------------------------------------------------
cat("STEP 3 - Intersection analysis at gene level", "\n")

## genes aligned but not overlapping known genes in annotation
wo_overlap_genes <- wo_overlap %>%
    group_by(seqnames_liftoff, strand_liftoff, query_gene_id, query_gene_biotype, query_gene_name) %>%
    summarise(max_tx_length = max(query_transcript_length),
           nb_tx = n()
    ) %>%
    distinct()

### determine query genes overlapping gtf
with_overlap_genes <- with_overlap %>%
    group_by(query_gene_id) %>% 
    summarise(nb_tx = n()) %>%
    distinct()
### add classification info: one-to-kwnown (at least one tx overlaps with gtf) or to-unknown gene (no tx overlapping with gtf)
wo_overlap_genes$classification <- ifelse(
    wo_overlap_genes$query_gene_id %in% with_overlap_genes$query_gene_id, 
    "one-to-known_gene", 
    "one-to-unknown_gene")

wo_overlap_genes <- wo_overlap_genes %>%
    subset(classification == "one-to-unknown_gene")

cat("Classification of query genes depending if at least one transcript overlaps target gtf", "\n")
print(table(wo_overlap_genes$classification, wo_overlap_genes$query_gene_biotype))

out_file_genes_wo <- paste0(species1,"_to_",species2,"_mapped_unknownGenes.txt")
write.table(x = wo_overlap_genes, file = out_file_genes_wo, quote = F, row.names = F, sep= "\t")
cat("File with",species1, "genes mapped on",species2, "but not overlapping with known genes in annotation: ", out_file_genes_wo, "\n")

## genes aligned and overlapping known genes in annotation

### add fraction overlap for each species
with_overlap$query_fraction_overlap <- with_overlap$overlap_length / with_overlap$query_transcript_length
with_overlap$target_fraction_overlap <- with_overlap$overlap_length / with_overlap$target_transcript_length

### group by gene and keep min/max value/fraction overlap obtained by tx
with_overlap_byGene <- with_overlap %>% 
    group_by(query_gene_id, query_gene_biotype, query_gene_name) %>%
    summarise(
        query_frac_overlap_min = min(query_fraction_overlap),
        query_frac_overlap_max = max(query_fraction_overlap),
        target_gene_id = paste(unique(target_gene_id), collapse = ","),
        target_gene_name = paste(unique(target_gene_name), collapse = ","),
        target_gene_biotype = paste(unique(target_gene_biotype), collapse = ","),
        target_frac_overlap_min = min(target_fraction_overlap),
        target_frac_overlap_max = max(target_fraction_overlap),
        strand_match = paste(unique(strand_match), collapse = ","),
        overlap_length_min = min(overlap_length),
        overlap_length_max = max(overlap_length))

with_overlap_byGene$strand_match <- ifelse(
    with_overlap_byGene$strand_match %in% c("antisense", "same_strand"),
    with_overlap_byGene$strand_match,
    "both"
)

### add query gene classification (nb target genes overlapping + classif)
with_overlap_byGene$nb_overlap_target_genes <- count.fields(textConnection(with_overlap_byGene$target_gene_id), 
                                                         sep = ",")

with_overlap_byGene <- with_overlap_byGene %>%
    mutate(query_gene_classif = case_when(nb_overlap_target_genes == 1 ~ "one_to_one",
                              nb_overlap_target_genes > 1 ~ "one_to_many",
                              TRUE ~ "other"))

cat("Classification of query genes overlapping target gtf", "\n")
print(table(with_overlap_byGene$query_gene_classif, with_overlap_byGene$query_gene_biotype))

out_file_genes_with <- paste0(species1,"_to_",species2,"_mapped_knownGenes.txt")
write.table(x = with_overlap_byGene, file = out_file_genes_with, quote = F, row.names = F, sep= "\t")
cat("File with ",species1, " genes mapped on ",species2, " and overlapping with known genes in annotation: ", out_file_genes_with, "\n")

