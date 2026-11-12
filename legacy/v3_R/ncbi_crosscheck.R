###############################################################################
##  NCBI RefSeq CROSS-CHECK  — catches coding genes Ensembl calls "lncRNA"
##  Your pipeline annotated overlaps against ENSEMBL only, which mislabels loci
##  that NCBI RefSeq annotates as protein-coding (e.g. SHROOM3, ADGRL3, ZNF391).
##  This re-checks every candidate against the NCBI RefSeq GFF too, and gives a
##  combined verdict: REJECT if coding in EITHER annotation.
###############################################################################
setwd("D:/dark_proteome_final")
suppressPackageStartupMessages({ library(openxlsx); library(rtracklayer); library(GenomicRanges) })
dir.create("ncbi", showWarnings=FALSE)
options(timeout=3600)

EVIDENCE <- "reports/dark_proteome_20aa_evidence.xlsx"

## ---- 1. download NCBI RefSeq GFF for ROS_Cfam_1.0 (GCF_014441545.1) --------
NCBI_GFF <- "ncbi/GCF_014441545.1_ROS_Cfam_1.0_genomic.gff"
if (!file.exists(NCBI_GFF)) {
  gz <- paste0(NCBI_GFF, ".gz")
  if (!file.exists(gz)) {
    cat("downloading NCBI RefSeq GFF ...\n")
    download.file("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/014/441/545/GCF_014441545.1_ROS_Cfam_1.0/GCF_014441545.1_ROS_Cfam_1.0_genomic.gff.gz",
                  gz, mode="wb", quiet=TRUE)
  }
  if (!requireNamespace("R.utils", quietly=TRUE)) install.packages("R.utils")
  R.utils::gunzip(gz, remove=FALSE, overwrite=TRUE)
}
stopifnot("NCBI GFF missing" = file.exists(NCBI_GFF))

## ---- 2. load candidate coordinates (from the evidence matrix) -------------
ev <- read.xlsx(EVIDENCE, sheet="Evidence_matrix")
beg <- if ("g_start" %in% names(ev)) ev$g_start else ev$start
end <- if ("g_end"   %in% names(ev)) ev$g_end   else ev$end
chr <- ev$chromosome; str <- ev$strand
ok  <- !is.na(chr) & !is.na(beg) & !is.na(end)
cat("candidates with coordinates:", sum(ok), "of", nrow(ev), "\n")

## ---- 3. CHROMOSOME NAME MAP: Ensembl ("1","2",..) -> NCBI ("NC_0518xx.1") --
## NCBI GFF uses RefSeq accessions as seqnames; your coords use Ensembl names.
## Build the map from the NCBI GFF's region features (they carry both).
gff <- rtracklayer::import(NCBI_GFF)
## region features list chromosome name in 'Name' or via 'chromosome' attr
reg <- gff[gff$type == "region"]
# Ensembl-style name is in the 'chromosome' or 'Name' mcol; accession is seqnames
acc  <- as.character(seqnames(reg))
# try to get the chromosome label
lab <- reg$Name
if (is.null(lab)) lab <- reg$chromosome
# keep only entries that have a simple chromosome label (1..38,X,MT)
map <- data.frame(ens = as.character(lab), ncbi = acc, stringsAsFactors=FALSE)
map <- map[!is.na(map$ens) & grepl("^[0-9XYMT]+$", map$ens), ]
map <- map[!duplicated(map$ens), ]
cat("chromosome map entries:", nrow(map), "(Ensembl name -> NCBI accession)\n")

## translate candidate chr (Ensembl) to NCBI accession
ncbi_chr <- map$ncbi[match(chr, map$ens)]

## ---- 4. NCBI gene features + biotype --------------------------------------
genes <- gff[gff$type == "gene"]
## NCBI stores biotype in 'gene_biotype'
gb <- genes$gene_biotype
gname <- genes$gene            # gene symbol
gr_genes <- GRanges(seqnames(genes), ranges(genes), strand="*")

## candidate GRanges in NCBI coordinates
cand <- GRanges(seqnames = ncbi_chr[ok],
                ranges   = IRanges(as.integer(beg[ok]), as.integer(end[ok])),
                strand   = "*")
hits <- findOverlaps(cand, gr_genes, ignore.strand=TRUE)

ncbi_gene <- rep(NA_character_, sum(ok))
ncbi_bt   <- rep("none",        sum(ok))
for (i in seq_along(cand)) {
  h <- subjectHits(hits)[queryHits(hits)==i]
  if (!length(h)) next
  bts <- gb[h]; nms <- gname[h]
  if (any(bts=="protein_coding", na.rm=TRUE)) {
    j <- which(bts=="protein_coding")[1]
    ncbi_gene[i] <- nms[j]; ncbi_bt[i] <- "protein_coding"
  } else {
    ncbi_gene[i] <- nms[1]; ncbi_bt[i] <- paste(unique(bts), collapse=";")
  }
}

## ---- 5. write NCBI columns back + COMBINED verdict ------------------------
ev$ncbi_gene    <- NA_character_; ev$ncbi_biotype <- NA_character_
ev$ncbi_gene[ok]    <- ncbi_gene
ev$ncbi_biotype[ok] <- ncbi_bt

ens_coding  <- grepl("protein_coding", ifelse(is.na(ev$overlap_biotype),"",ev$overlap_biotype))
ncbi_coding <- grepl("protein_coding", ifelse(is.na(ev$ncbi_biotype),"",ev$ncbi_biotype))

ev$coding_either <- ens_coding | ncbi_coding
ev$final_verdict <- ifelse(is.na(ev$peptide_novel) | ev$peptide_novel==FALSE, "KNOWN (peptide)",
                           ifelse(ev$coding_either, "REJECT (coding in Ensembl or NCBI)",
                                  "CANDIDATE (non-coding both)"))

cat("\n=== COMBINED (Ensembl + NCBI) VERDICT ===\n")
print(table(ev$final_verdict))

## clean candidate shortlist
cols <- intersect(c("orf_id","aa_length","n_tissues","max_runs","peptide_novel",
                    "overlap_biotype","ncbi_biotype","ncbi_gene","conserved","final_verdict"), names(ev))
cand_tbl <- ev[ev$final_verdict=="CANDIDATE (non-coding both)", cols]
cand_tbl <- cand_tbl[order(-cand_tbl$n_tissues, -cand_tbl$max_runs, -cand_tbl$aa_length), ]
cat("\n=== CLEAN CANDIDATES (non-coding in BOTH annotations) ===\n")
print(cand_tbl, row.names=FALSE)

## also show what NCBI rescued (Ensembl said lncRNA, NCBI says coding)
rescued <- ev[!ens_coding & ncbi_coding & ev$peptide_novel==TRUE,
              c("orf_id","overlap_biotype","ncbi_gene","ncbi_biotype")]
cat("\n=== CAUGHT BY NCBI (Ensembl=lncRNA but NCBI=coding -> now rejected) ===\n")
print(rescued, row.names=FALSE)

wb <- loadWorkbook(EVIDENCE)
if ("Evidence_matrix" %in% names(wb)) removeWorksheet(wb, "Evidence_matrix")
addWorksheet(wb,"Evidence_matrix"); writeData(wb,"Evidence_matrix",ev)
saveWorkbook(wb, EVIDENCE, overwrite=TRUE)
write.xlsx(cand_tbl, "reports/CLEAN_candidates.xlsx", overwrite=TRUE)
cat("\nUpdated evidence matrix + wrote reports/CLEAN_candidates.xlsx\n")