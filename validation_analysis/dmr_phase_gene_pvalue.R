#!/usr/bin/env Rscript

library(ggplot2)
library(tidyverse)
library(data.table)
library(grid)
library(ggvenn)
library(janitor)
library(glue)
library(argparse)

######################################################################################################
#### argument parsing
parser <- ArgumentParser(description="a R script to do analyze the differentially methylation data.")

parser$add_argument("-d", "--dmr", type="character", dest="dmr_list", 
                    help="provide a list with complete file names separated by comma", required=F)

parser$add_argument("-m", "--methyl", type="character", dest="methyl_list", 
                    help="provide a list with complete file names separated by comma", required=F)

parser$add_argument("-p", "--path", type="character", dest="path", 
                    help="provide the absolute path for input files", required=TRUE)

parser$add_argument("-a", "--agroup", type="character", dest="a_group", 
                    help="provide a group name", required=TRUE)

parser$add_argument("-b", "--bgroup", type="character", dest="b_group", 
                    help="provide b group name", required=TRUE)

parser$add_argument("-c", "--condition", type="character", dest="condition", 
                    help="provide b group name", required=TRUE)

parser$add_argument("-s", "--set", type="character", dest="data_set", 
                    help="provide data set type", required=TRUE)

parser$add_argument("-o", "--outdir", type="character", dest="out_dir", 
                    help="provide the absolute path for output files", required=TRUE)

args <- parser$parse_args()

path <- args$path
dmr_list <- unlist(strsplit(args$dmr_list, ","))
methyl_list <- unlist(strsplit(args$methyl_list, ","))
a_group <- args$a_group
b_group <- args$b_group
con <- args$condition
dataset <- args$data_set
out_dir <- args$out_dir

######################################################################################################
############ function blocks
#### Function to read dmr files
read_dmr_files <- function(pattern, dmr_file_list) {
  files <- dmr_list[grepl(pattern, dmr_file_list)]
  data <- lapply(files, function(f) {
    df <- fread(f, sep = "\t", header = TRUE)
    return(df)
  })
  rbindlist(data)
}

#### Function to read methyl files
read_methyl_files <- function(pattern, methyl_file_list, colname) {
  files <- methyl_list[grepl(pattern, methyl_file_list)]
  data <- lapply(files, function(f) {
    df <- fread(f, sep = "\t", header = FALSE)
    setnames(df, colname)
    return(df)
  })
  rbindlist(data)
}

position_clean <- function(pileup_position_filter, dmr_tsv_file) {
  position_file <- pileup_position_filter %>%
    mutate(specific_id_strand = paste(pileup_position_filter$chrom,pileup_position_filter$end, pileup_position_filter$strand, sep = "_"),
           specific_id = paste(pileup_position_filter$chrom, pileup_position_filter$end, sep = "_")) %>%
    select(specific_id_strand, specific_id, strand)
  
  position_file_tidy <- position_file[!duplicated(position_file$specific_id_strand), ]
  position_file_tidy_vec <- sub("_[+-]$", "", position_file_tidy$specific_id_strand)
  
  duplicate_id <- position_file_tidy_vec[duplicated(position_file_tidy_vec)]
  print(duplicate_id)
  position_file_tidy_clean <- position_file_tidy %>%
    filter(!(specific_id %in% duplicate_id))
  
  dmr_tsv_file <- dmr_tsv_file %>%
    mutate(specific_id = paste(dmr_tsv_file$chrom, dmr_tsv_file$end, sep = "_"))
  
  dmr_tsv_file_clean <- dmr_tsv_file %>%
    filter(!(specific_id %in% duplicate_id))
  
  dmr_tsv_file_clean <- dmr_tsv_file_clean %>%
    select(-strand) %>%
    inner_join(position_file_tidy_clean, by = join_by(specific_id), keep = F, unmatched = "drop")
  
  unique_end_num <- length(unique(dmr_tsv_file_clean$specific_id))
  dmr_tsv_file_clean <- dmr_tsv_file_clean %>%
    mutate(differential=ifelse(balanced_map_pvalue <= 0.05 & abs(cohen_h) >= 0.2, "yes", "no"),
           group=paste(a_group, b_group, sep="_"),
           condition=con)
  
  return(dmr_tsv_file_clean)
}


#### a function to change the comparison group name in the dmr file
group_name_change <- function(dmr_tsv_file_clean, a_group_name, b_group_name) {
  names(dmr_tsv_file_clean) <- gsub("^a_", paste0(a_group_name, "_", sep=""), colnames(dmr_tsv_file_clean))
  names(dmr_tsv_file_clean) <- gsub("^b_", paste0(b_group_name, "_", sep=""), colnames(dmr_tsv_file_clean))
  
  names(dmr_tsv_file_clean) <- gsub("_a_", paste0("_", a_group_name, "_", sep=""), colnames(dmr_tsv_file_clean))
  names(dmr_tsv_file_clean) <- gsub("_b_", paste0("_", b_group_name, "_", sep=""), colnames(dmr_tsv_file_clean))
  
  return(dmr_tsv_file_clean)
}

#### functions to combine methylated sites and genes
deferential_methylated_form <- function(dmr_methyl_position, gff_tsv){
  differential_methyl_from <- as.data.frame(matrix(ncol = length(c(names(dmr_methyl_position), "id", "locus_tag", "symbol")), nrow = 0))
  
  colnames(differential_methyl_from) <- c(names(dmr_methyl_position), "id", "locus_tag", "symbol")
  
  for (i in seq_len(nrow(gff_tsv))){
    
    seqname <- as.character(gff_tsv[i,"seqname"])
    gff_start <- as.integer(gff_tsv[i,"start"])
    gff_end <- as.integer(gff_tsv[i,"end"])
    gff_strand <- as.character(gff_tsv[i,"strand"])
    gff_id <- as.character(gff_tsv[i,"id"])
    gff_locus <- as.character(gff_tsv[i,"locus_tag"])
    gff_symbol <- as.character(gff_tsv[i,"symbol"])
    
    inter_form <- subset(dmr_methyl_position,
                         dmr_methyl_position$chrom == seqname & dmr_methyl_position$end >= gff_start & dmr_methyl_position$strand == gff_strand &
                           dmr_methyl_position$end <= gff_end)
    
    inter_form$id <- gff_id
    inter_form$locus_tag <- gff_locus
    inter_form$symbol <- gff_symbol
    
    differential_methyl_from <- rbind(differential_methyl_from, inter_form)
    
  }
  
  differential_methyl_from <- differential_methyl_from %>%
    mutate(condition=con)
  
  return(differential_methyl_from)
}

##### a functions to extract the lowest p-valur for genes with differential 
methyl_form_differential_sort_pavlue <- function(differential_methyl_from, pileup_position, specified_region){
  differential_methyl_from_pvalue <- differential_methyl_from %>%
    group_by(chrom, id, condition) %>%
    slice_min(., balanced_map_pvalue) %>%
    group_by(chrom, condition) %>%
    distinct(id, .keep_all = T) %>%
    select(chrom, start, end, specific_id_strand, balanced_map_pvalue, cohen_h, differential, id, locus_tag, symbol, condition)
  
  pileup_position <- pileup_position %>%
    mutate(specific_id_strand = paste(chrom, end, strand, sep = "_"))
  
  differential_methyl_from_pvalue <- differential_methyl_from_pvalue[differential_methyl_from_pvalue$specific_id_strand %in% pileup_position$specific_id_strand,]
  
  pileup_unique <- unique(pileup_position[, c("specific_id_strand", "modified_base_code")])
  
  differential_methyl_from_pvalue_clean <- merge(differential_methyl_from_pvalue, pileup_unique[, c("specific_id_strand","modified_base_code")], by = "specific_id_strand", all.x = TRUE) %>%
    mutate(region = specified_region)
  
  return(differential_methyl_from_pvalue_clean)
}

#### a function to clean up the modified code in the piple up file
position_clean_modified_code <- function(pileup_position, raw_modified_code, clean_modified_code) {
  pileup_position_clean_modified_code <- pileup_position %>%
    filter(modified_base_code == raw_modified_code) %>%
    mutate(modified_base_code = clean_modified_code)
  return(pileup_position_clean_modified_code)
}

######################################################################################################
#### loading file
cat("\n now is loading file...🏃 ")
setwd(path)

colname_methyl <- c("chrom",	"start",	"end",	"modified_base_code",
                    "score",	"strand",	"valid_coverage",	"modified_number",	
                    "phase",	"replicate",	"condition")

m6a_dmr_files <- read_dmr_files("_m6a_", dmr_list)
m4c_dmr_files <- read_dmr_files("_m4c_", dmr_list)
m5c_dmr_files <- read_dmr_files("_m5c_", dmr_list)

m6a_metyhl_files <- read_methyl_files("m6a", methyl_list, colname_methyl)
m4c_metyhl_files <- read_methyl_files("m4c", methyl_list, colname_methyl)
m5c_metyhl_files <- read_methyl_files("m5c", methyl_list, colname_methyl)

ecoli_mRNA_promoter <- fread("GCF_000008865.2_ASM886v2_genomic_gene_extended_promoter.tsv", sep="\t", header = F)
ecoli_gene_body <- fread("GCF_000008865.2_ASM886v2_genomic_gene.tsv", sep="\t", header = F)
ecoli_nonregulate <- fread("GCF_000008865.2_ASM886v2_genomic_nonregulatory.bed", sep="\t", header = F)

colnames(ecoli_mRNA_promoter) <- c("seqname", "feature", "start", "end", "strand", "id", "locus_tag", "symbol")
colnames(ecoli_gene_body) <- c("seqname", "feature", "start", "end", "strand", "id", "locus_tag", "symbol")
colnames(ecoli_nonregulate) <- c("seqname", "start", "end")
options(scipen=999)

#######################################################################################################################################
#### perform file cleaning and methylation region tageting only on differential sties
cat("\n now is cleaning file and generating new differential files...🧹 \n")
ecoli_m6A_position_clean <- position_clean(m6a_metyhl_files, m6a_dmr_files)
ecoli_m4C_position_clean <- position_clean(m4c_metyhl_files, m4c_dmr_files)
ecoli_m5C_position_clean <- position_clean(m5c_metyhl_files, m5c_dmr_files)

ecoli_m6A_position_clean <- group_name_change(ecoli_m6A_position_clean, a_group, b_group)
ecoli_m4C_position_clean <- group_name_change(ecoli_m4C_position_clean, a_group, b_group)
ecoli_m5C_position_clean <- group_name_change(ecoli_m5C_position_clean, a_group, b_group)

ecoli_m4C_position_clean_differntial <- ecoli_m4C_position_clean
ecoli_m5C_position_clean_differntial <- ecoli_m5C_position_clean
ecoli_m6A_position_clean_differntial <- ecoli_m6A_position_clean

comparison <- paste(a_group, b_group, sep="_")


###### m6a
promoter_methyl_from_differential_m6A <- deferential_methylated_form(ecoli_m6A_position_clean_differntial, ecoli_mRNA_promoter)
gene_body_methyl_from_differential_m6A <- deferential_methylated_form(ecoli_m6A_position_clean_differntial, ecoli_gene_body)


####### m4c
promoter_methyl_from_differential_m4C <- deferential_methylated_form(ecoli_m4C_position_clean_differntial, ecoli_mRNA_promoter)
gene_body_methyl_from_differential_m4C <- deferential_methylated_form(ecoli_m4C_position_clean_differntial, ecoli_gene_body)


###### m5c
promoter_methyl_from_differential_m5C <- deferential_methylated_form(ecoli_m5C_position_clean_differntial, ecoli_mRNA_promoter)
gene_body_methyl_from_differential_m5C <- deferential_methylated_form(ecoli_m5C_position_clean_differntial, ecoli_gene_body)



#######################################################################################################################################
#### identify gene with differenetial methylation and only extract the lowest p-value loci
###### m6a
promoter_sort_pavlue_m6A <- methyl_form_differential_sort_pavlue(promoter_methyl_from_differential_m6A, 
                                                                 position_clean_modified_code(m6a_metyhl_files, "a", "m6A"), "promoter")
gene_body_sort_pavlue_m6A <- methyl_form_differential_sort_pavlue(gene_body_methyl_from_differential_m6A, 
                                                                  position_clean_modified_code(m6a_metyhl_files, "a", "m6A"), "gene_body")

###### m4c
promoter_sort_pavlue_m4C <- methyl_form_differential_sort_pavlue(promoter_methyl_from_differential_m4C, 
                                                                 position_clean_modified_code(m4c_metyhl_files, "21839", "m4C"), "promoter")
gene_body_sort_pavlue_m4C <- methyl_form_differential_sort_pavlue(gene_body_methyl_from_differential_m4C, 
                                                                  position_clean_modified_code(m4c_metyhl_files, "21839", "m4C"), "gene_body")

###### m5c
promoter_sort_pavlue_m5C <- methyl_form_differential_sort_pavlue(promoter_methyl_from_differential_m5C, 
                                                                 position_clean_modified_code(m5c_metyhl_files, "m", "m5C"), "promoter")
gene_body_sort_pavlue_m5C <- methyl_form_differential_sort_pavlue(gene_body_methyl_from_differential_m5C, 
                                                                  position_clean_modified_code(m5c_metyhl_files, "m", "m5C"), "gene_body")


gene_body_sort_pavluel_combine <- rbind(gene_body_sort_pavlue_m6A, gene_body_sort_pavlue_m4C, gene_body_sort_pavlue_m5C)

promoter_sort_pavluel_combine <- rbind(promoter_sort_pavlue_m6A, promoter_sort_pavlue_m4C, promoter_sort_pavlue_m5C)

gene_body_promoter_sort_pavluel_combine <- rbind(gene_body_sort_pavluel_combine,
                                                 promoter_sort_pavluel_combine) %>%
  mutate(set=dataset)

write_tsv(gene_body_promoter_sort_pavluel_combine, 
          file=glue("{out_dir}/{a_group}_{b_group}_{con}_{dataset}_gene_body_promoter_sort_pvalue.tsv"))
