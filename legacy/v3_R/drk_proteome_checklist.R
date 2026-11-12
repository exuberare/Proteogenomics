###############################################################################
##  DARK PROTEOME — MASTER CHECKLIST BUILDER
##  Appends every candidate (with BLAST + tBLASTn detail) to one running
##  Excel tracker. Re-run any time: new candidates are ADDED, your manual
##  Remark/Decision edits on existing rows are PRESERVED.
##
##  One row per ORF with: PRIDE id, ORF id, gene, transcript, locus, length,
##  tissue detection, novelty verdict, genomic context, and a blank Remark +
##  Decision column for you to fill in by hand.
###############################################################################

## ======================= CONFIG =========================
PROJECT    <- "D:/dark_proteome_final"
META_CSV   <- "orfs/lncRNA_orf_metadata_20aa.csv"
EVIDENCE   <- "reports/dark_proteome_20aa_evidence.xlsx"   # from the analysis script
CHECKLIST  <- "reports/DARK_PROTEOME_CHECKLIST.xlsx"       # the running tracker
PRIDE_ID   <- "PXD052932"                                  # dataset accession
ASSEMBLY   <- "Ensembl ROS_Cfam_1.0 (GCF_014441545.1)"
## ========================================================

setwd(PROJECT)
suppressPackageStartupMessages(library(openxlsx))
stopifnot("Run the analysis script first (evidence file missing)." = file.exists(EVIDENCE))

ev   <- read.xlsx(EVIDENCE, sheet="Evidence_matrix")
meta <- read.csv(META_CSV, stringsAsFactors=FALSE)

## ---- assemble the checklist rows (one per ORF) ----------------------------
runcols <- grep("_runs$", names(ev), value=TRUE)
tissue_of <- function(r){                      # which tissues detected this ORF
  t <- sub("_runs$","", runcols[as.numeric(r[runcols]) > 0]); paste(t, collapse=", ") }
new <- data.frame(
  PRIDE_ID        = PRIDE_ID,
  ORF_ID          = ev$orf_id,
  gene_id         = ev$gene_id,
  gene_name       = meta$gene_name[match(ev$orf_id, meta$orf_id)],
  transcript_id   = ev$transcript_id,
  aa_length       = ev$aa_length,
  num_exons       = ev$num_exons,
  tissues         = apply(ev, 1, tissue_of),
  n_tissues       = ev$n_tissues,
  max_runs        = ev$max_runs,
  reproducible    = ev$max_runs >= 2,
  peptide_novel   = ev$peptide_novel,
  blast_verdict   = ev$verdict,
  blast_pident    = if("blast_pident" %in% names(ev)) ev$blast_pident else NA,
  blast_evalue    = if("blast_evalue" %in% names(ev)) ev$blast_evalue else NA,
  chromosome      = ev$chromosome,
  start           = if("g_start" %in% names(ev)) ev$g_start else NA,
  end             = if("g_end"   %in% names(ev)) ev$g_end   else NA,
  strand          = ev$strand,
  locus           = ev$locus,
  overlaps_known_CDS = ev$overlaps_known_CDS,
  overlapping_gene   = if("overlapping_gene" %in% names(ev)) ev$overlapping_gene else NA,
  assembly        = ASSEMBLY,
  aa_sequence     = meta$aa_sequence[match(ev$orf_id, meta$orf_id)],
  date_added      = as.character(Sys.Date()),
  ## ---- manual columns (you fill these in) ----
  Repeat_region   = "",          # yes/no — from genome browser (RepeatMasker)
  Overlaps_coding_gene_browser = "",  # yes/no — what the browser actually shows
  Remark          = "",          # your free-text note
  Decision        = ""           # CONSIDER / REJECT / MAYBE
)

## ---- merge with existing checklist (preserve manual edits) ----------------
MANUAL <- c("Repeat_region","Overlaps_coding_gene_browser","Remark","Decision")
if (file.exists(CHECKLIST)) {
  old <- read.xlsx(CHECKLIST, sheet=1)
  # carry over your manual entries for ORFs already in the list
  for (mc in MANUAL) if (mc %in% names(old)) {
    idx <- match(new$ORF_ID, old$ORF_ID)
    carried <- old[[mc]][idx]
    new[[mc]] <- ifelse(!is.na(carried) & carried != "", carried, new[[mc]])
  }
  # keep any old ORFs that aren't in the new run (nothing lost)
  extra <- old[!(old$ORF_ID %in% new$ORF_ID), , drop=FALSE]
  if (nrow(extra)) {
    for (cn in setdiff(names(new), names(extra))) extra[[cn]] <- NA
    new <- rbind(new, extra[, names(new)])
  }
  n_added <- sum(!(new$ORF_ID %in% old$ORF_ID))
  cat("existing checklist:", nrow(old), "rows | new ORFs added:", n_added, "\n")
} else {
  cat("creating new checklist:", nrow(new), "ORFs\n")
}

## order: strongest first (novel + reproducible on top)
new$peptide_novel[is.na(new$peptide_novel)] <- FALSE
new <- new[order(-as.integer(new$peptide_novel), -new$max_runs, -new$aa_length), ]

## ---- write with formatting ------------------------------------------------
wb <- createWorkbook()
addWorksheet(wb, "Checklist")
writeData(wb, "Checklist", new, withFilter=TRUE)
freezePane(wb, "Checklist", firstActiveRow=2, firstActiveCol=3)
# header style
hs <- createStyle(fgFill="#1F3864", fontColour="#FFFFFF", textDecoration="bold", halign="center")
addStyle(wb, "Checklist", hs, rows=1, cols=1:ncol(new), gridExpand=TRUE)
# highlight the manual columns so they stand out
mc_idx <- which(names(new) %in% MANUAL)
ms <- createStyle(fgFill="#FFF2CC")
addStyle(wb, "Checklist", ms, rows=2:(nrow(new)+1), cols=mc_idx, gridExpand=TRUE, stack=TRUE)
setColWidths(wb, "Checklist", cols=1:ncol(new), widths="auto")
setColWidths(wb, "Checklist", cols=which(names(new)=="aa_sequence"), widths=40)
saveWorkbook(wb, CHECKLIST, overwrite=TRUE)

cat("\nwritten", CHECKLIST, "—", nrow(new), "ORFs total\n")
cat("Manual columns to fill (highlighted yellow): ", paste(MANUAL, collapse=", "), "\n")