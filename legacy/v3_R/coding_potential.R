###############################################################################
##  CODING-POTENTIAL CHECK for candidate dark ORFs
##  For ORFs with NO conservation (can't be validated by homology), this asks:
##  does the sequence look like genuine protein-coding sequence?
##
##  Produces inputs for CPAT / CPC2 (web tools) and computes quick in-house
##  metrics (Kozak context, ORF/transcript coverage, GC, aa composition flags).
###############################################################################
setwd("D:/dark_proteome_final")
suppressPackageStartupMessages({ library(Biostrings); library(openxlsx) })
dir.create("coding_potential", showWarnings=FALSE)

META_CSV <- "orfs/lncRNA_orf_metadata_20aa.csv"
EVIDENCE <- "reports/dark_proteome_20aa_evidence.xlsx"
TX_FA    <- "transcripts/lncRNA_transcripts.fa"   # source transcripts (nucleotide)

ev   <- read.xlsx(EVIDENCE, sheet="Evidence_matrix")
meta <- read.csv(META_CSV, stringsAsFactors=FALSE)

## ---- choose which ORFs to test -------------------------------------------
## focus on your real candidates: reproducible OR flagged CONSIDER/novel.
## edit this filter as you like; default = all candidates in the evidence matrix
ids <- unique(ev$orf_id)
sub <- meta[match(ids, meta$orf_id), ]
cat("ORFs to assess:", nrow(sub), "\n")

## ---- 1. export NUCLEOTIDE ORF sequences for CPAT/CPC2 ---------------------
## CPAT/CPC2 want the ORF nucleotide (CDS) sequence. Reconstruct it from the
## transcript using nt_start/nt_end in the metadata.
tx <- readDNAStringSet(TX_FA); names(tx) <- sub("\\s.*","",names(tx))
getnt <- function(i){
  t <- sub$transcript_id[i]
  if (!t %in% names(tx)) return(NA_character_)
  s <- as.character(tx[[t]])
  en <- min(nchar(s), sub$nt_end[i])          # guard transcript end
  substr(s, sub$nt_start[i], en)
}
nt <- vapply(seq_len(nrow(sub)), getnt, character(1))
names(nt) <- sub$orf_id
nt <- nt[!is.na(nt) & nchar(nt) >= 18]
orf_nt <- DNAStringSet(nt)
writeXStringSet(orf_nt, "coding_potential/candidate_ORFs_nucleotide.fasta")
cat("wrote nucleotide ORFs for CPAT/CPC2:", length(orf_nt),
    "-> coding_potential/candidate_ORFs_nucleotide.fasta\n")

## also protein FASTA (for CPC2, which can take either / for records)
aa <- AAStringSet(setNames(sub$aa_sequence, sub$orf_id))
writeXStringSet(aa, "coding_potential/candidate_ORFs_protein.fasta")

## ---- 2. quick in-house coding-support metrics -----------------------------
## (a) Kozak context around the ATG: strong if -3 is A/G and +4 is G
kozak <- function(i){
  t <- sub$transcript_id[i]; if (!t %in% names(tx)) return(NA_character_)
  s <- as.character(tx[[t]]); p <- sub$nt_start[i]
  m3 <- if (p-3 >= 1) substr(s,p-3,p-3) else "N"   # -3 position
  p4 <- if (p+3 <= nchar(s)) substr(s,p+3,p+3) else "N" # +4 position
  strong <- (m3 %in% c("A","G")) & (p4 == "G")
  paste0("-3=",m3," +4=",p4, if (strong) " (strong)" else " (weak)")
}
## (b) ORF coverage of its transcript (longer coverage ~ more coding-like)
cov <- round(100 * (sub$nt_end - sub$nt_start) /
               width(tx)[match(sub$transcript_id, names(tx))], 1)
## (c) GC content of the ORF nucleotide
gcpct <- round(100 * letterFrequency(orf_nt, "GC", as.prob=TRUE)[match(sub$orf_id, names(orf_nt))], 1)
## (d) low-complexity / homopolymer flag on the protein
lowcomp <- vapply(sub$aa_sequence, function(a){
  f <- table(strsplit(a,"")[[1]]); max(f)/nchar(a) > 0.4 }, logical(1))

qc <- data.frame(orf_id=sub$orf_id, aa_length=sub$aa_length,
                 kozak=vapply(seq_len(nrow(sub)), kozak, character(1)),
                 orf_cover_pct=cov, gc_pct=as.numeric(gcpct),
                 low_complexity=lowcomp,
                 conserved=ev$conserved[match(sub$orf_id, ev$orf_id)],
                 stringsAsFactors=FALSE)
qc <- qc[order(qc$conserved, -qc$aa_length), ]
print(qc, row.names=FALSE)
write.xlsx(qc, "coding_potential/coding_support_metrics.xlsx", overwrite=TRUE)
cat("\nwrote coding_potential/coding_support_metrics.xlsx\n")

cat("\n=== NEXT: run CPAT or CPC2 (web) on the nucleotide FASTA ===\n")
cat("  CPC2  : http://cpc2.gao-lab.org/  (upload candidate_ORFs_nucleotide.fasta)\n")
cat("  CPAT  : https://wlcb.oit.uic.edu/cpat/ (choose the closest model: Human/Mouse/Fly/Zebrafish)\n")
cat("  -> CPC2 labels each as 'coding' or 'noncoding' with a coding probability.\n")
cat("  -> Download the CPC2/CPAT result table and we'll merge it into the evidence matrix.\n")