=========================================================
  PANDA Demo Dataset Overview
=========================================================
This dataset provides simulated bisulfite sequencing data (NGS & Sanger)
to demonstrate the capabilities of PANDA at single-base resolution.

[1] Reference Sequence (Upload this first!)
---------------------------------------------------------
File: NGS_Ref_Targets.fasta (or Sanger_Ref_Targets.fasta)
 * mm10_H19_DMR: Mouse H19 ICR for Exp A & B.
 * mm10_Nanog_Promoter: Real genomic region (Nanog locus) for Exp C.

=========================================================
[Experiment Details & How to use in PANDA]
=========================================================

[Exp A] Allele-Specific Methylation (ASM) & Loss of Imprinting (LOI)
---------------------------------------------------------
Target: mm10_H19_DMR
 * Normal_Rep: Demonstrates normal imprinting. Expect a bimodal distribution.
 * LOI_Rep: Demonstrates a disease state where both alleles are aberrantly hypermethylated.

[Exp B] Haplotype Phasing (SNP-based allele separation)
---------------------------------------------------------
Target: mm10_H19_DMR
The file 'Hetero_Phased_Rep' contains a 50:50 mixture of two haplotypes (A and B).
Initially, PANDA will show a bimodal histogram. To demonstrate phasing, enter the
following 21-bp Motif to cleanly filter out Haplotype B and isolate Haplotype A:

 -> Motif A (Isolates Unmethylated allele): TGTATAAATGATTGATTTTTT
 -> Motif B (Isolates Methylated allele)  : TGTATAAATGGTTGATTTTTT

[Exp C] Group Comparison (Nanog TF Footprint Remodeling)
---------------------------------------------------------
Target: mm10_Nanog_Promoter
 * WT vs KO: Demonstrates dynamic epigenetic remodeling.
PANDA's 'Difference Plot' will beautifully reveal a mix of up-regulated (red) and
down-regulated (blue) CpG sites within the exact same region.
