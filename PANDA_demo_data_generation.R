############################################################
# PANDA Demo Data Generator
############################################################

set.seed(11)

library(BSgenome.Mmusculus.UCSC.mm10)
library(Biostrings)

ROOT_DIR <- "PANDA_Demo"
GT_DIR   <- file.path(ROOT_DIR, "GroundTruth") 

DIRS <- list(
  SANGER_REF   = file.path(ROOT_DIR, "Sanger", "Reference"),
  SANGER_MULTI = file.path(ROOT_DIR, "Sanger", "Multi"),
  SANGER_SINGLE= file.path(ROOT_DIR, "Sanger", "Single"),
  NGS_REF      = file.path(ROOT_DIR, "NGS", "Reference"),
  NGS_FASTQ    = file.path(ROOT_DIR, "NGS", "FASTQ"),
  NGS_FASTA    = file.path(ROOT_DIR, "NGS", "FASTA")
)

if (dir.exists(ROOT_DIR)) unlink(ROOT_DIR, recursive=TRUE)
for (d in DIRS) dir.create(d, recursive=TRUE)
dir.create(GT_DIR, recursive=TRUE)

message("Extracting mm10 reference sequences...")
mm10 <- BSgenome.Mmusculus.UCSC.mm10

get_dense_amplicon <- function(genome, chr, start_pos, end_pos, width = 350) {
  wide_seq <- getSeq(genome, chr, start_pos, end_pos)
  cpg_hits <- start(matchPattern("CG", wide_seq))
  cpg_counts <- sapply(1:(length(wide_seq)-width), function(i) sum(cpg_hits >= i & cpg_hits < i+width))
  best_idx <- which.max(cpg_counts)
  best_seq <- subseq(wide_seq, best_idx, best_idx+width-1)
  message(sprintf(">>> SUCCESS: Found dense region on %s with %d CpGs (width: %dbp) <<<", chr, max(cpg_counts), width))
  return(list(seq = best_seq, count = max(cpg_counts)))
}

message("Scanning H19 locus for a highly dense CpG amplicon...")
h19_res <- get_dense_amplicon(mm10, "chr7", 142570000, 142585000, 350)
SEQ_H19 <- h19_res$seq
n_cpg_h19 <- h19_res$count

message("Scanning the Nanog promoter locus for a highly dense CpG amplicon...")
nanog_res <- get_dense_amplicon(mm10, "chr6", 122700000, 122715000, 350)
SEQ_NANOG <- nanog_res$seq
n_cpg_target <- nanog_res$count

writeXStringSet(DNAStringSet(list(mm10_H19_DMR = SEQ_H19, mm10_Nanog_Promoter = SEQ_NANOG)), file.path(DIRS$NGS_REF, "NGS_Ref_Targets.fasta"))
writeXStringSet(DNAStringSet(list(mm10_H19_DMR = SEQ_H19, mm10_Nanog_Promoter = SEQ_NANOG)), file.path(DIRS$SANGER_REF, "Sanger_Ref_Targets.fasta"))

# Exp B Setup: Dynamically introduce a central SNP (A/G)
s_ref_h19 <- strsplit(as.character(SEQ_H19), "")[[1]]
cpg_idx_h19 <- which(s_ref_h19 == "C" & c(s_ref_h19[-1], "") == "G")
non_cpg_c_idx <- setdiff(which(s_ref_h19 == "C"), cpg_idx_h19)

center_pos <- length(s_ref_h19) / 2
target_snp_pos <- non_cpg_c_idx[which.min(abs(non_cpg_c_idx - center_pos))]

SEQ_A <- replaceAt(SEQ_H19, IRanges(target_snp_pos, target_snp_pos), DNAString("A"))
SEQ_B <- replaceAt(SEQ_H19, IRanges(target_snp_pos, target_snp_pos), DNAString("G"))

motif_start <- max(1, target_snp_pos - 10)
motif_end <- min(length(s_ref_h19), target_snp_pos + 10)

s_ref_A <- strsplit(as.character(SEQ_A), "")[[1]]
s_ref_A[setdiff(which(s_ref_A == "C"), cpg_idx_h19)] <- "T"
MOTIF_A <- paste(s_ref_A[motif_start:motif_end], collapse="")

s_ref_B <- strsplit(as.character(SEQ_B), "")[[1]]
s_ref_B[setdiff(which(s_ref_B == "C"), cpg_idx_h19)] <- "T"
MOTIF_B <- paste(s_ref_B[motif_start:motif_end], collapse="")

simulate_reads <- function(base_seq, mean_meth, num_reads, err_rate) {
  s_ref <- strsplit(as.character(base_seq), "")[[1]]
  seq_len <- length(s_ref)
  cpg_indices <- which(s_ref == "C" & c(s_ref[-1], "") == "G")
  non_cpg_c_indices <- setdiff(which(s_ref == "C"), cpg_indices)
  
  if(length(mean_meth) == 1) {
    mean_meth <- rep(mean_meth, length(cpg_indices))
  } 
  
  reads <- character(num_reads)
  actual_meths <- numeric(num_reads)
  bases <- c("A", "T", "G", "C")
  
  for(i in 1:num_reads) {
    s_read <- s_ref
    s_read[non_cpg_c_indices] <- "T"
    
    read_meth_sites <- numeric(length(cpg_indices))
    for(j in seq_along(cpg_indices)) {
      alpha_param <- max(mean_meth[j] * 50, 0.5)
      beta_param  <- max((1 - mean_meth[j]) * 50, 0.5)
      read_meth_sites[j] <- rbeta(1, alpha_param, beta_param)
    }
    actual_meths[i] <- mean(read_meth_sites)
    
    meth_probs <- runif(length(cpg_indices))
    convert_to_t <- cpg_indices[meth_probs > read_meth_sites]
    s_read[convert_to_t] <- "T"
    
    if(err_rate > 0) {
      err_idx <- which(runif(seq_len) < err_rate)
      for(idx in err_idx) s_read[idx] <- sample(bases[bases != s_read[idx]], 1)
    }
    reads[i] <- paste(s_read, collapse="")
  }
  return(list(reads = reads, meths = actual_meths))
}

simulate_ngs <- function(base_seq, mean_meth_vec, props, sample_name, experiment, replicate, folder) {
  coverage <- 3000
  seq_len <- length(base_seq)
  props <- props / sum(props)
  covs <- round(coverage * props)
  if(sum(covs) != coverage) covs[1] <- covs[1] + (coverage - sum(covs))
  
  all_reads <- character(0); all_meths <- numeric(0); all_alleles <- character(0)
  
  for(k in 1:length(props)) {
    if(covs[k] > 0) {
      sub_meth <- if(is.list(mean_meth_vec)) mean_meth_vec[[k]] else mean_meth_vec[k]
      sim <- simulate_reads(base_seq, sub_meth, covs[k], err_rate = 0.005)
      all_reads <- c(all_reads, sim$reads)
      all_meths <- c(all_meths, sim$meths)
      all_alleles <- c(all_alleles, rep(paste0("Subpopulation_", k), covs[k]))
    }
  }
  
  idx <- sample(1:coverage)
  all_reads <- all_reads[idx]; all_meths <- all_meths[idx]; all_alleles <- all_alleles[idx]
  
  gt_df <- data.frame(Sample = sample_name, ReadID = paste0("Read_", 1:coverage), sampled_methylation = all_meths, allele = all_alleles)
  write.csv(gt_df, file.path(GT_DIR, paste0(sample_name, "_GT.csv")), row.names=FALSE)
  
  fq_dir <- file.path(DIRS$NGS_FASTQ, folder)
  dir.create(fq_dir, recursive=TRUE, showWarnings=FALSE)
  fq_path <- file.path(fq_dir, paste0(sample_name, ".fastq.gz"))
  lines <- character(coverage * 4)
  lines[c(TRUE, FALSE, FALSE, FALSE)] <- paste0("@Read_", 1:coverage)
  lines[c(FALSE, TRUE, FALSE, FALSE)] <- all_reads
  lines[c(FALSE, FALSE, TRUE, FALSE)] <- "+"
  lines[c(FALSE, FALSE, FALSE, TRUE)] <- rep(paste(rep("I", seq_len), collapse=""), coverage)
  
  con <- gzfile(fq_path, "w")
  writeLines(lines, con)
  close(con)
  
  fa_dir <- file.path(DIRS$NGS_FASTA, folder)
  dir.create(fa_dir, recursive=TRUE, showWarnings=FALSE)
  fa_path <- file.path(fa_dir, paste0(sample_name, ".fasta"))
  fa_lines <- character(coverage * 2)
  fa_lines[c(TRUE, FALSE)] <- paste0(">Read_", 1:coverage)
  fa_lines[c(FALSE, TRUE)] <- all_reads
  writeLines(fa_lines, fa_path)
}

simulate_mixed_ngs <- function(seq_A, seq_B, mean_A, mean_B, sample_name, folder) {
  coverage <- 3000
  cov_A <- coverage / 2
  cov_B <- coverage / 2
  
  sim_A <- simulate_reads(seq_A, mean_A, cov_A, err_rate = 0.005)
  sim_B <- simulate_reads(seq_B, mean_B, cov_B, err_rate = 0.005)
  
  all_reads <- c(sim_A$reads, sim_B$reads)
  all_meths <- c(sim_A$meths, sim_B$meths)
  all_alleles <- c(rep("Haplotype_A_Unmeth", cov_A), rep("Haplotype_B_Meth", cov_B))
  
  idx <- sample(1:coverage)
  all_reads <- all_reads[idx]; all_meths <- all_meths[idx]; all_alleles <- all_alleles[idx]
  
  gt_df <- data.frame(Sample = sample_name, ReadID = paste0("Read_", 1:coverage), sampled_methylation = all_meths, allele = all_alleles)
  write.csv(gt_df, file.path(GT_DIR, paste0(sample_name, "_GT.csv")), row.names=FALSE)
  
  seq_len <- length(seq_A)
  fq_dir <- file.path(DIRS$NGS_FASTQ, folder)
  dir.create(fq_dir, recursive=TRUE, showWarnings=FALSE)
  fq_path <- file.path(fq_dir, paste0(sample_name, ".fastq.gz"))
  lines <- character(coverage * 4)
  lines[c(TRUE, FALSE, FALSE, FALSE)] <- paste0("@Read_", 1:coverage)
  lines[c(FALSE, TRUE, FALSE, FALSE)] <- all_reads
  lines[c(FALSE, FALSE, TRUE, FALSE)] <- "+"
  lines[c(FALSE, FALSE, FALSE, TRUE)] <- rep(paste(rep("I", seq_len), collapse=""), coverage)
  
  con <- gzfile(fq_path, "w")
  writeLines(lines, con)
  close(con)
}

simulate_sanger <- function(base_seq, mean_meth_vec, props, sample_name, experiment, replicate, folder, single=FALSE) {
  clones <- 16
  props <- props / sum(props)
  covs <- round(clones * props)
  if(sum(covs) != clones) covs[1] <- covs[1] + (clones - sum(covs))
  
  all_reads <- character(0); all_meths <- numeric(0)
  for(k in 1:length(props)) {
    if(covs[k] > 0) {
      sub_meth <- if(is.list(mean_meth_vec)) mean_meth_vec[[k]] else mean_meth_vec[k]
      sim <- simulate_reads(base_seq, sub_meth, covs[k], err_rate = 0.001)
      all_reads <- c(all_reads, sim$reads); all_meths <- c(all_meths, sim$meths)
    }
  }
  
  idx <- sample(1:clones)
  all_reads <- all_reads[idx]; all_meths <- all_meths[idx]
  fasta_names <- paste0("Clone_", 1:clones)
  
  gt_df <- data.frame(Sample = sample_name, clone = fasta_names, sampled_methylation = all_meths)
  write.csv(gt_df, file.path(GT_DIR, paste0(sample_name, "_GT.csv")), row.names=FALSE)
  
  out_dir <- if(single) DIRS$SANGER_SINGLE else file.path(DIRS$SANGER_MULTI, folder)
  dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
  names(all_reads) <- fasta_names
  writeXStringSet(DNAStringSet(all_reads), file.path(out_dir, paste0(sample_name, ".fasta")))
}

simulate_mixed_sanger <- function(seq_A, seq_B, mean_A, mean_B, sample_name, folder) {
  clones <- 16
  cov_A <- clones / 2
  cov_B <- clones / 2
  
  sim_A <- simulate_reads(seq_A, mean_A, cov_A, err_rate = 0.001)
  sim_B <- simulate_reads(seq_B, mean_B, cov_B, err_rate = 0.001)
  
  all_reads <- c(sim_A$reads, sim_B$reads)
  all_meths <- c(sim_A$meths, sim_B$meths)
  
  idx <- sample(1:clones)
  all_reads <- all_reads[idx]; all_meths <- all_meths[idx]
  fasta_names <- paste0("Clone_", 1:clones)
  
  gt_df <- data.frame(Sample = sample_name, clone = fasta_names, sampled_methylation = all_meths)
  write.csv(gt_df, file.path(GT_DIR, paste0(sample_name, "_GT.csv")), row.names=FALSE)
  
  out_dir <- file.path(DIRS$SANGER_MULTI, folder)
  dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
  names(all_reads) <- fasta_names
  writeXStringSet(DNAStringSet(all_reads), file.path(out_dir, paste0(sample_name, ".fasta")))
}

message("Generating Demo Data (Exp A, B, C) and Ground Truth...")

wt_meth_pattern <- rep(0.98, n_cpg_target)
ko_meth_pattern <- rep(0.98, n_cpg_target)

if (n_cpg_target >= 3) {
  mid_start <- floor(n_cpg_target * 0.3) + 1
  mid_end   <- ceiling(n_cpg_target * 0.7)
  wt_meth_pattern[mid_start:mid_end] <- 0.02
  ko_meth_pattern[1:(mid_start - 1)] <- 0.02 
}

for(i in 1:3){
  # Exp A
  simulate_ngs(SEQ_H19, c(0.02, 0.98), c(0.5, 0.5), paste0("NGS_H19_Normal_Rep",i), "ExpA_ASM", i, "ExpA_ASM")
  simulate_ngs(SEQ_H19, c(0.98), c(1.0), paste0("NGS_H19_LOI_Rep",i), "ExpA_ASM", i, "ExpA_ASM")
  simulate_sanger(SEQ_H19, c(0.02, 0.98), c(0.5, 0.5), paste0("Sanger_H19_Normal_Rep",i), "ExpA_ASM", i, "ExpA_ASM")
  simulate_sanger(SEQ_H19, c(0.98), c(1.0), paste0("Sanger_H19_LOI_Rep",i), "ExpA_ASM", i, "ExpA_ASM")
  
  simulate_mixed_ngs(SEQ_A, SEQ_B, c(0.02), c(0.98), paste0("NGS_H19_Hetero_Phased_Rep",i), "ExpB_Phasing")
  simulate_mixed_sanger(SEQ_A, SEQ_B, c(0.02), c(0.98), paste0("Sanger_H19_Hetero_Phased_Rep",i), "ExpB_Phasing")
  
  # Exp C
  simulate_ngs(SEQ_NANOG, list(wt_meth_pattern), c(1.0), paste0("NGS_Nanog_WT_Rep",i), "ExpC_Group", i, "ExpC_Group")
  simulate_sanger(SEQ_NANOG, list(wt_meth_pattern), c(1.0), paste0("Sanger_Nanog_WT_Rep",i), "ExpC_Group", i, "ExpC_Group")
  simulate_ngs(SEQ_NANOG, list(ko_meth_pattern), c(1.0), paste0("NGS_Nanog_KO_Rep",i), "ExpC_Group", i, "ExpC_Group")
  simulate_sanger(SEQ_NANOG, list(ko_meth_pattern), c(1.0), paste0("Sanger_Nanog_KO_Rep",i), "ExpC_Group", i, "ExpC_Group")
}

simulate_sanger(SEQ_H19, c(0.02, 0.98), c(0.5, 0.5), "Sanger_Single_H19_ASM", "Single_Test", 1, "", single=TRUE)

readme_text <- c(
  "=========================================================",
  "  PANDA Demo Dataset Overview",
  "=========================================================",
  "This dataset provides simulated bisulfite sequencing data (NGS & Sanger)",
  "to demonstrate the capabilities of PANDA at single-base resolution.",
  "",
  "[1] Reference Sequence (Upload this first!)",
  "---------------------------------------------------------",
  "File: NGS_Ref_Targets.fasta (or Sanger_Ref_Targets.fasta)",
  " * mm10_H19_DMR: Mouse H19 ICR for Exp A & B.",
  " * mm10_Nanog_Promoter: Real genomic region (Nanog locus) for Exp C.",
  "",
  "=========================================================",
  "[Experiment Details & How to use in PANDA]",
  "=========================================================",
  "",
  "[Exp A] Allele-Specific Methylation (ASM) & Loss of Imprinting (LOI)",
  "---------------------------------------------------------",
  "Target: mm10_H19_DMR",
  " * Normal_Rep: Demonstrates normal imprinting. Expect a bimodal distribution.",
  " * LOI_Rep: Demonstrates a disease state where both alleles are aberrantly hypermethylated.",
  "",
  "[Exp B] Haplotype Phasing (SNP-based allele separation)",
  "---------------------------------------------------------",
  "Target: mm10_H19_DMR",
  "The file 'Hetero_Phased_Rep' contains a 50:50 mixture of two haplotypes (A and B).",
  "Initially, PANDA will show a bimodal histogram. To demonstrate phasing, enter the",
  "following 21-bp Motif to cleanly filter out Haplotype B and isolate Haplotype A:",
  "",
  paste0(" -> Motif A (Isolates Unmethylated allele): ", MOTIF_A),
  paste0(" -> Motif B (Isolates Methylated allele)  : ", MOTIF_B),
  "",
  "[Exp C] Group Comparison (Nanog TF Footprint Remodeling)",
  "---------------------------------------------------------",
  "Target: mm10_Nanog_Promoter",
  " * WT vs KO: Demonstrates dynamic epigenetic remodeling.",
  "PANDA's 'Difference Plot' will beautifully reveal a mix of up-regulated (red) and",
  "down-regulated (blue) CpG sites within the exact same region."
)

writeLines(readme_text, file.path(ROOT_DIR, "README.txt"))
message("\n>>> SUCCESS: Demo Data, GroundTruth CSVs, and Detailed README generated! <<<")