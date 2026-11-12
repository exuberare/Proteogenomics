###############################################################################
##  AUTO-ANNOTATE ORF GENOMIC OVERLAP  (replaces manual browser tooltip checks)
##  For every candidate ORF (using its tBLASTn genomic coordinates), reports:
##    - overlapping gene(s) + biotype  (protein_coding vs lncRNA/non-coding)
##    - whether the ORF is EXONIC or INTRONIC in the overlapping transcript(s)
##    - an auto-decision following your manual rule
##  Then updates the evidence matrix and checklist with these columns.
###############################################################################
setwd("D:/dark_proteome_final")
suppressPackageStartupMessages({
  library(openxlsx); library(rtracklayer); library(GenomicRanges)
})

GTF_PATH <- "reference_genome/genome/Canis_lupus_familiaris.ROS_Cfam_1.0.115.gtf"
EVIDENCE <- "reports/dark_proteome_20aa_evidence.xlsx"

## ---- load candidates with their genomic coordinates ----------------------
ev <- read.xlsx(EVIDENCE, sheet="Evidence_matrix")
beg <- if ("g_start" %in% names(ev)) ev$g_start else ev$start
end <- if ("g_end"   %in% names(ev)) ev$g_end   else ev$end
chr <- ev$chromosome
ok  <- !is.na(chr) & !is.na(beg) & !is.na(end)
cat("ORFs with coordinates:", sum(ok), "of", nrow(ev), "\n")

orf_gr <- GRanges(seqnames = chr[ok],
                  ranges   = IRanges(as.integer(beg[ok]), as.integer(end[ok])),
                  strand   = ifelse(ev$strand[ok] %in% c("+","-"), ev$strand[ok], "*"))
names(orf_gr) <- ev$orf_id[ok]

## ---- load GTF: genes (with biotype) and exons ----------------------------
gtf <- rtracklayer::import(GTF_PATH)
# normalise seqnames of GTF to match your ORF coords (both should be "1","2",...,"X")
# (your tBLASTn used the Ensembl genome, so chr names already match the GTF)
genes <- gtf[gtf$type == "gene"]
exons <- gtf[gtf$type == "exon"]

## ---- 1. overlap with GENES (any strand) -> biotype -----------------------
## we take ALL overlapping genes, then decide by the "worst" biotype
## (protein_coding overlap dominates -> reject)
hits <- findOverlaps(orf_gr, genes, ignore.strand = TRUE)
gene_biotype_of <- function(i){
  h <- subjectHits(hits)[queryHits(hits) == i]
  if (!length(h)) return(c(gene="", biotype="none"))
  bt <- genes$gene_biotype[h]
  nm <- genes$gene_name[h]; nm[is.na(nm)] <- genes$gene_id[h][is.na(nm)]
  # prioritise protein_coding if present
  if (any(bt == "protein_coding")) {
    j <- which(bt == "protein_coding")[1]
    return(c(gene=nm[j], biotype="protein_coding"))
  }
  c(gene=nm[1], biotype=paste(unique(bt), collapse=";"))
}

## ---- 2. EXONIC vs INTRONIC ------------------------------------------------
## exonic if the ORF overlaps any exon of an overlapping gene; else (overlaps a
## gene but no exon) = intronic; else = intergenic/none
exon_hit <- overlapsAny(orf_gr, exons, ignore.strand = TRUE)
gene_hit <- overlapsAny(orf_gr, genes, ignore.strand = TRUE)
overlap_type <- ifelse(exon_hit, "exonic",
                       ifelse(gene_hit, "intronic", "intergenic"))

## ---- assemble per-ORF annotation -----------------------------------------
n <- length(orf_gr)
gb <- t(vapply(seq_len(n), gene_biotype_of, c(gene="", biotype="")))
ann <- data.frame(orf_id = names(orf_gr),
                  overlapping_gene = gb[,"gene"],
                  overlap_biotype  = gb[,"biotype"],
                  overlap_type     = overlap_type,
                  stringsAsFactors = FALSE)

## ---- 3. AUTO-DECISION following your manual rule --------------------------
auto_decision <- function(bt, ot){
  if (bt == "protein_coding")                 return("REJECT (coding gene)")
  if (ot == "intergenic")                     return("CONSIDER (no gene overlap)")
  if (grepl("lncRNA|ncRNA|processed|TEC|antisense", bt) && ot == "exonic")
    return("CONSIDER (non-coding, exonic)")
  if (ot == "intronic")                       return("MAYBE (intronic)")
  return("MAYBE (review)")
}
ann$auto_decision <- mapply(auto_decision, ann$overlap_biotype, ann$overlap_type)

cat("\n=== AUTO-CLASSIFICATION ===\n")
print(table(ann$auto_decision))

## ---- merge back into evidence matrix + write ------------------------------
m <- match(ev$orf_id, ann$orf_id)
ev$overlapping_gene_auto <- ann$overlapping_gene[m]
ev$overlap_biotype       <- ann$overlap_biotype[m]
ev$overlap_type          <- ann$overlap_type[m]
ev$auto_decision         <- ann$auto_decision[m]

wb <- loadWorkbook(EVIDENCE)
if ("Evidence_matrix" %in% names(wb)) removeWorksheet(wb, "Evidence_matrix")
addWorksheet(wb, "Evidence_matrix"); writeData(wb, "Evidence_matrix", ev)
saveWorkbook(wb, EVIDENCE, overwrite = TRUE)

## also print the candidate view
show <- ev[order(ev$auto_decision, -ev$max_runs), 
           c("orf_id","aa_length","max_runs","peptide_novel",
             "overlapping_gene_auto","overlap_biotype","overlap_type","auto_decision")]
print(show, row.names = FALSE)
cat("\nUpdated", EVIDENCE, "with auto-overlap columns.\n")
cat("Re-run your checklist script to pull these into the checklist.\n")