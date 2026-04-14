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

parser$add_argument("-g", "--gp", type="character", dest="growth_phase", 
                    help="provide b group name", required=TRUE)

parser$add_argument("-o", "--outdir", type="character", dest="out_dir", 
                    help="provide the absolute path for output files", required=TRUE)

args <- parser$parse_args()

path <- args$path
dmr_list <- unlist(strsplit(args$dmr_list, ","))
methyl_list <- unlist(strsplit(args$methyl_list, ","))
a_group <- args$a_group
b_group <- args$b_group
growth_phase <- args$growth_phase
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
           phase=growth_phase)
  
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
    
    #ifelse(inter_form$specific_id != "", print(paste(inter_form$specific_id, "is done", sep=" ")), "")

  }
  
  differential_methyl_from <- differential_methyl_from %>%
    mutate(group=paste(a_group, b_group, sep="_"))
  
  return(differential_methyl_from)
}


deferential_methylated_form_none_regulatory <- function(dmr_methyl_position, gff_tsv){
  differential_methyl_from <- as.data.frame(matrix(ncol = length(c(names(dmr_methyl_position))), nrow = 0))

  colnames(differential_methyl_from) <- c(names(dmr_methyl_position))
  
  for (i in seq_len(nrow(gff_tsv))){
    
    seqname <- as.character(gff_tsv[i,"seqname"])
    gff_start <- as.integer(gff_tsv[i,"start"])
    gff_end <- as.integer(gff_tsv[i,"end"])
    
    inter_form <- subset(dmr_methyl_position,
                         dmr_methyl_position$chrom == seqname & dmr_methyl_position$end > gff_start & dmr_methyl_position$end < gff_end)
    
    differential_methyl_from <- rbind(differential_methyl_from, inter_form)
    
    
  }
  
  differential_methyl_from <- differential_methyl_from %>%
    mutate(group=paste(a_group, b_group, sep="_"))
  
  return(differential_methyl_from)
}

#### a function to clean up the differential form with methylated sites and genes
#### clean up differential form
differential_form_clean_up <- function(differential_form) {
  differential_form_clean <- differential_form %>%
    select(chrom, end, exp_pct_modified, sta_pct_modified, symbol)
  
  differential_form_clean <- differential_form_clean %>%
    rename(exponential = exp_pct_modified,
           stationary = sta_pct_modified)
  
  differential_form_clean_longer <- differential_form_clean %>%
    pivot_longer(cols = c("exponential", "stationary"), names_to = "phase", values_to = "modified_percentage")
  
  differential_form_clean_summary <- differential_form_clean_longer %>%
    group_by(chrom, phase, symbol) %>%
    summarize(mean_modified_percentage = mean(modified_percentage), 
              sd_modified_percentage = sd(modified_percentage), .groups= "drop")
  return(differential_form_clean_summary)
}

#### functions to identify differential modified position at each replicate with gene information
differntial_position_replicate <- function(differential_methyl_from, pileup_position, specified_region){
  differ_position_list <- differential_methyl_from %>%
    select(specific_id_strand, id, locus_tag, symbol)
  
  
  pileup_position <- pileup_position %>%
    mutate(specific_id_strand = paste(chrom, end, strand, sep = "_"))
  
  
  differential_methyl_region_position_clean <- pileup_position[pileup_position$specific_id_strand %in% differ_position_list$specific_id_strand,]
  
  
  differential_methyl_region_position <- merge(differ_position_list, differential_methyl_region_position_clean, by = "specific_id_strand", all.x = TRUE)
  

  
  differential_methyl_region_position_summary <- differential_methyl_region_position %>%
    group_by(chrom, start, end, strand, modified_base_code, phase, condition, id, locus_tag, symbol) %>%
    summarise(mean_modified_fraction = mean(100*modified_number/valid_coverage), 
              sd_modified_fraction = sd(100*modified_number/valid_coverage), 
              region = specified_region, .groups = "drop")
  
  differential_methyl_region_position_summary[, c("start", "end")] <- lapply(differential_methyl_region_position_summary[, c("start", "end")], as.character)
  
  return(differential_methyl_region_position_summary)
}

differntial_position_replicate_none_regulatory <- function(differential_methyl_from, pileup_position, specified_region){
  differ_position_list <- differential_methyl_from %>%
    select(specific_id_strand)
  
  pileup_position <- pileup_position %>%
    mutate(specific_id_strand = paste(chrom, end, strand, sep = "_"))
  
  differential_methyl_region_position_clean <- pileup_position[pileup_position$specific_id_strand %in% differ_position_list$specific_id_strand,]
  
  differential_methyl_region_position <- merge(differ_position_list, differential_methyl_region_position_clean, by = "specific_id_strand", all.x = TRUE)
  
  differential_methyl_region_position_summary <- differential_methyl_region_position %>%
    group_by(chrom, start, end, strand, modified_base_code, phase, condition) %>%
    summarise(mean_modified_fraction = mean(100*modified_number/valid_coverage), 
              sd_modified_fraction = sd(100*modified_number/valid_coverage), 
              region = specified_region, .groups = "drop")
  
  differential_methyl_region_position_summary[, c("start", "end")] <- lapply(differential_methyl_region_position_summary[, c("start", "end")], as.character)
  
  return(differential_methyl_region_position_summary)
}

#### a function to clean up the modified code in the pile up file
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

m6a_metyhl_files <- read_methyl_files("m6a", methyl_list, colname_methyl) %>%
  filter(phase==growth_phase)
m4c_metyhl_files <- read_methyl_files("m4c", methyl_list, colname_methyl) %>%
  filter(phase==growth_phase)
m5c_metyhl_files <- read_methyl_files("m5c", methyl_list, colname_methyl) %>%
  filter(phase==growth_phase)

ecoli_mRNA_promoter <- fread("GCF_000008865.2_ASM886v2_genomic_gene_extended_promoter.tsv", sep="\t", header = F)
ecoli_gene_body <- fread("GCF_000008865.2_ASM886v2_genomic_gene.tsv", sep="\t", header = F)
ecoli_nonregulate <- fread("GCF_000008865.2_ASM886v2_genomic_nonregulatory.bed", sep="\t", header = F)

colnames(ecoli_mRNA_promoter) <- c("seqname", "feature", "start", "end", "strand", "id", "locus_tag", "symbol")
colnames(ecoli_gene_body) <- c("seqname", "feature", "start", "end", "strand", "id", "locus_tag", "symbol")
colnames(ecoli_nonregulate) <- c("seqname", "start", "end")
options(scipen=999)


if (growth_phase == "exponential") {
  phase_short <- "exp"
} else if (growth_phase == "stationary") {
  phase_short <- "sta"
} else {
  cat("\nUnidentified phase term\n")
}

#######################################################################################################################################
#### perform file cleaning and methylation region tageting only on differential sties
cat("\n now is cleaning file and generating new differential files...🧹 \n")
ecoli_m6A_position_clean <- position_clean(m6a_metyhl_files, m6a_dmr_files)
ecoli_m4C_position_clean <- position_clean(m4c_metyhl_files, m4c_dmr_files)
ecoli_m5C_position_clean <- position_clean(m5c_metyhl_files, m5c_dmr_files)

ecoli_m6A_position_clean <- group_name_change(ecoli_m6A_position_clean, a_group, b_group)
ecoli_m4C_position_clean <- group_name_change(ecoli_m4C_position_clean, a_group, b_group)
ecoli_m5C_position_clean <- group_name_change(ecoli_m5C_position_clean, a_group, b_group)

ecoli_m4C_position_clean_differntial <- ecoli_m4C_position_clean[ecoli_m4C_position_clean$differential=="yes",]
ecoli_m5C_position_clean_differntial <- ecoli_m5C_position_clean[ecoli_m5C_position_clean$differential=="yes",]
ecoli_m6A_position_clean_differntial <- ecoli_m6A_position_clean[ecoli_m6A_position_clean$differential=="yes",]

comparison <- paste(a_group, b_group, sep="_")

####### m4c
promoter_methyl_from_differential_m4C <- deferential_methylated_form(ecoli_m4C_position_clean_differntial, ecoli_mRNA_promoter)
gene_body_methyl_from_differential_m4C <- deferential_methylated_form(ecoli_m4C_position_clean_differntial, ecoli_gene_body)
neither_methyl_from_differential_m4C <- deferential_methylated_form_none_regulatory(ecoli_m4C_position_clean_differntial, ecoli_nonregulate)

###### m5c
promoter_methyl_from_differential_m5C <- deferential_methylated_form(ecoli_m5C_position_clean_differntial, ecoli_mRNA_promoter)
gene_body_methyl_from_differential_m5C <- deferential_methylated_form(ecoli_m5C_position_clean_differntial, ecoli_gene_body)
neither_methyl_from_differential_m5C <- deferential_methylated_form_none_regulatory(ecoli_m5C_position_clean_differntial, ecoli_nonregulate)

###### m6a
promoter_methyl_from_differential_m6A <- deferential_methylated_form(ecoli_m6A_position_clean_differntial, ecoli_mRNA_promoter)
gene_body_methyl_from_differential_m6A <- deferential_methylated_form(ecoli_m6A_position_clean_differntial, ecoli_gene_body)
neither_methyl_from_differential_m6A <- deferential_methylated_form_none_regulatory(ecoli_m6A_position_clean_differntial, ecoli_nonregulate)

ecoli_position_clean_differntial <- rbind(ecoli_m4C_position_clean_differntial, ecoli_m5C_position_clean_differntial, ecoli_m6A_position_clean_differntial)
write_tsv(ecoli_position_clean_differntial, 
          file=glue("{out_dir}/{a_group}_{b_group}_{phase_short}_position_clean_differntial.tsv"))

ecoli_position_clean_differntial_nonregulatory <- rbind(neither_methyl_from_differential_m4C, neither_methyl_from_differential_m5C, neither_methyl_from_differential_m6A)
write_tsv(ecoli_position_clean_differntial_nonregulatory, 
          file=glue("{out_dir}/{a_group}_{b_group}_{phase_short}_position_clean_differntial_nonregulatory.tsv"))


#######################################################################################################################################
#### identify the differential methylated positions around regularoty and non-regulatory regions
#### m6a
cat("\n now is identifying the methylation positions...🔍 ")
promoter_m6a_position_differential <- differntial_position_replicate(promoter_methyl_from_differential_m6A,
                                                                     position_clean_modified_code(m6a_metyhl_files, "a", "m6A"), "promoter")

gene_body_m6a_position_differential <- differntial_position_replicate(gene_body_methyl_from_differential_m6A,
                                                                      position_clean_modified_code(m6a_metyhl_files, "a", "m6A"), "gene_body")

none_regulatory_m6a_position_differential <- differntial_position_replicate_none_regulatory(neither_methyl_from_differential_m6A,
                                                                                            position_clean_modified_code(m6a_metyhl_files, "a", "m6A"), "neither")

##### m4c
promoter_m4c_position_differential <- differntial_position_replicate(promoter_methyl_from_differential_m4C,
                                                                     position_clean_modified_code(m4c_metyhl_files, "21839", "m4C"),
                                                                     "promoter")

gene_body_m4c_position_differential <- differntial_position_replicate(gene_body_methyl_from_differential_m4C,
                                                                      position_clean_modified_code(m4c_metyhl_files, "21839", "m4C"), "gene_body")

none_regulatory_m4c_position_differential <- differntial_position_replicate_none_regulatory(gene_body_methyl_from_differential_m4C,
                                                                                            position_clean_modified_code(m4c_metyhl_files, "21839", "m4C"), "neither")

#### m5c
promoter_m5c_position_differential <- differntial_position_replicate(promoter_methyl_from_differential_m5C,
                                                                     position_clean_modified_code(m5c_metyhl_files, "m", "m5C"), "promoter")
gene_body_m5c_position_differential <- differntial_position_replicate(gene_body_methyl_from_differential_m5C,
                                                                      position_clean_modified_code(m5c_metyhl_files, "m", "m5C"), "gene_body")

none_regulatory_m5c_position_differential <- differntial_position_replicate_none_regulatory(gene_body_methyl_from_differential_m5C,
                                                                                            position_clean_modified_code(m5c_metyhl_files, "m", "m5C"), "neither")

#### file combination
gene_body_position_differntial_combine <- gene_body_m6a_position_differential %>%
  rbind(gene_body_m4c_position_differential, gene_body_m5c_position_differential)

promoter_position_differntial_combine <- promoter_m6a_position_differential %>%
  rbind(promoter_m4c_position_differential, promoter_m5c_position_differential)

gene_body_promoter_position_differntial_combine <- rbind(gene_body_position_differntial_combine,
                                                         promoter_position_differntial_combine)

write_tsv(gene_body_promoter_position_differntial_combine, 
          file=glue("{out_dir}/{a_group}_{b_group}_{phase_short}_diffferential_metyhl_position_gene_summary.tsv"))

none_regulatory_position_differntial_combine <- rbind(none_regulatory_m6a_position_differential, none_regulatory_m4c_position_differential, none_regulatory_m5c_position_differential)

write_tsv(none_regulatory_position_differntial_combine, 
          file=glue("{out_dir}/{a_group}_{b_group}_{phase_short}_diffferential_metyhl_none_regulatory_position_summary.tsv"))

#######################################################################################################################################
#### print the gene id of m6a, m4c, m5c
### m6A
# promoter
cat("\n now is generating the gene id...💻 ")
gene_id_promoter_m6A <- sub("^GeneID:", "", unique(promoter_methyl_from_differential_m6A$id))
gene_symbol_promoter_m6A <- unique(promoter_methyl_from_differential_m6A$symbol)

write_lines(gene_id_promoter_m6A, glue("{out_dir}/{a_group}_{b_group}_{phase_short}_gene_list_promoter_m6A.tsv"))
write_lines(gene_symbol_promoter_m6A, glue("{out_dir}/{a_group}_{b_group}_{phase_short}_gene_symbol_list_promoter_m6A.tsv"))

# gene body
gene_id_gene_body_m6A <- sub("^GeneID:", "", unique(gene_body_methyl_from_differential_m6A$id))
gene_symbol_gene_body_m6A <- unique(gene_body_methyl_from_differential_m6A$symbol)

write_lines(gene_id_gene_body_m6A, glue("{out_dir}/{a_group}_{b_group}_{phase_short}_gene_list_gene_body_m6A.tsv"))
write_lines(gene_symbol_gene_body_m6A, glue("{out_dir}/{a_group}_{b_group}_{phase_short}_gene_symbol_list_gene_body_m6A.tsv"))

################################################################################
### m4C
# promoter
gene_id_promoter_m4C <- sub("^GeneID:", "", unique(promoter_methyl_from_differential_m4C$id))
gene_symbol_promoter_m4C <- unique(promoter_methyl_from_differential_m4C$symbol)

write_lines(gene_id_promoter_m4C, glue("{out_dir}/{a_group}_{b_group}_{phase_short}_gene_list_promoter_m4C.tsv"))
write_lines(gene_symbol_promoter_m4C, glue("{out_dir}/{a_group}_{b_group}_{phase_short}_gene_symbol_list_promoter_m4C.tsv"))

# gene body
gene_id_gene_body_m4C <- sub("^GeneID:", "", unique(gene_body_methyl_from_differential_m4C$id))
gene_symbol_gene_body_m4C <- unique(gene_body_methyl_from_differential_m4C$symbol)

write_lines(gene_id_gene_body_m4C, glue("{out_dir}/{a_group}_{b_group}_{phase_short}_gene_list_gene_body_m4C.tsv"))
write_lines(gene_symbol_gene_body_m4C, glue("{out_dir}/{a_group}_{b_group}_{phase_short}_gene_symbol_list_gene_body_m4C.tsv"))

################################################################################
### m5C
# promoter
gene_id_promoter_m5C <- sub("^GeneID:", "", unique(promoter_methyl_from_differential_m5C$id))
gene_symbol_promoter_m5C <- unique(promoter_methyl_from_differential_m5C$symbol)

write_lines(gene_id_promoter_m5C, glue("{out_dir}/{a_group}_{b_group}_{phase_short}_gene_list_promoter_m5C.tsv"))
write_lines(gene_symbol_promoter_m5C, glue("{out_dir}/{a_group}_{b_group}_{phase_short}_gene_symbol_list_promoter_m5C.tsv"))

# gene body
gene_id_gene_body_m5C <- sub("^GeneID:", "", unique(gene_body_methyl_from_differential_m5C$id))
gene_symbol_gene_body_m5C <- unique(gene_body_methyl_from_differential_m5C$symbol)

write_lines(gene_id_gene_body_m5C, glue("{out_dir}/{a_group}_{b_group}_{phase_short}_gene_list_gene_body_m5C.tsv"))
write_lines(gene_symbol_gene_body_m5C, glue("{out_dir}/{a_group}_{b_group}_{phase_short}_gene_symbol_list_gene_body_m5C.tsv"))


################################################################################
### all methylation forms
# promoter
gene_id_promoter_all_methyl <- unique(c(gene_id_promoter_m4C, gene_id_promoter_m5C, gene_id_promoter_m6A))
gene_symbol_promoter_all_methyl <- unique(c(gene_symbol_promoter_m4C, gene_symbol_promoter_m5C, gene_symbol_promoter_m6A))

gene_id_gene_body_all_methyl <- unique(c(gene_id_gene_body_m4C, gene_id_gene_body_m5C, gene_id_gene_body_m6A))
gene_symbol_gene_body_all_methyl <- unique(c(gene_symbol_gene_body_m4C, gene_symbol_gene_body_m5C, gene_symbol_gene_body_m6A))

write_lines(gene_id_promoter_all_methyl, glue("{out_dir}/{a_group}_{b_group}_{phase_short}_gene_list_promoter_all_methyl.tsv"))
write_lines(gene_symbol_promoter_all_methyl, glue("{out_dir}/{a_group}_{b_group}_{phase_short}_gene_symbol_list_promoter_all_methyl.tsv"))

write_lines(gene_id_gene_body_all_methyl, glue("{out_dir}/{a_group}_{b_group}_{phase_short}_gene_list_gene_body_all_methyl.tsv"))
write_lines(gene_symbol_gene_body_all_methyl, glue("{out_dir}/{a_group}_{b_group}_{phase_short}_gene_symbol_list_gene_body_all_methyl.tsv"))

cat("\n the analysis of this round is done! 🍻\n")
cat("\n")
