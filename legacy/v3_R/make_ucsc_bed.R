###############################################################################
##  EXPORT CANDIDATE ORFs AS A UCSC CUSTOM TRACK (BED)
##  Turns the tBLASTn genomic coordinates in the evidence matrix into a BED
##  file you upload to the UCSC browser (ROS_Cfam_1.0 / GCF_014441545.1 hub)
##  to see ALL candidate ORFs as one track alongside RefSeq / RepeatMasker.
###############################################################################
setwd("D:/dark_proteome_final")
suppressPackageStartupMessages(library(openxlsx))

EVIDENCE <- "reports/dark_proteome_20aa_evidence.xlsx"
OUT_BED  <- "reports/candidate_ORFs_ucsc.bed"

ev <- read.xlsx(EVIDENCE, sheet="Evidence_matrix")

## coordinate columns (handle either naming)
chr <- ev$chromosome
beg <- if ("g_start" %in% names(ev)) ev$g_start else ev$start
end <- if ("g_end"   %in% names(ev)) ev$g_end   else ev$end
str <- ev$strand
ok  <- !is.na(chr) & !is.na(beg) & !is.na(end)
cat("candidates with genomic coordinates:", sum(ok), "of", nrow(ev), "\n")
ev <- ev[ok,]; chr<-chr[ok]; beg<-beg[ok]; end<-end[ok]; str<-str[ok]

## ---- chromosome name handling --------------------------------------------
## Your tBLASTn ran against the Ensembl genome, so seqnames are Ensembl style
## ("1","2",...,"X","MT") OR RefSeq accessions ("NC_051805.1"). The UCSC
## ROS_Cfam_1.0 hub uses "chr1","chr2",... So prefix bare numbers/letters with
## "chr". RefSeq accessions (NC_0518xx) are already what the hub shows as chrN
## via its alias — if your chr looks like an accession, see note at bottom.
ucsc_chr <- ifelse(grepl("^(chr)", chr), chr,
                   ifelse(grepl("^[0-9XYM]+$|^MT$", chr),
                          paste0("chr", sub("^MT$","M",chr)), chr))

## ---- strand / colour by decision-ish status -------------------------------
strand <- ifelse(str %in% c("+","-"), str, ".")
# colour: novel = green, known/other = grey (purely visual)
novel <- !is.na(ev$peptide_novel) & ev$peptide_novel == TRUE
col   <- ifelse(novel, "0,158,115", "150,150,150")   # green / grey

## ---- BED12-ish (BED6 + thick + rgb) ---------------------------------------
bed <- data.frame(
  chrom      = ucsc_chr,
  chromStart = as.integer(beg) - 1L,     # BED is 0-based start
  chromEnd   = as.integer(end),
  name       = ev$orf_id,
  score      = pmin(1000, ev$max_runs * 150),  # darker = more reproducible
  strand     = strand,
  thickStart = as.integer(beg) - 1L,
  thickEnd   = as.integer(end),
  itemRgb    = col,
  stringsAsFactors = FALSE)
bed <- bed[order(bed$chrom, bed$chromStart),]

## ---- write with the UCSC track header -------------------------------------
hdr <- 'track name="DarkORFs" description="Candidate lncRNA dark ORFs" visibility=pack itemRgb="On"'
writeLines(hdr, OUT_BED)
write.table(bed, OUT_BED, sep="\t", quote=FALSE, row.names=FALSE,
            col.names=FALSE, append=TRUE)

cat("written", OUT_BED, "with", nrow(bed), "ORFs\n")
cat("chromosomes used:", paste(sort(unique(bed$chrom)), collapse=", "), "\n")
cat("\nUpload: UCSC genome-euro -> My Data -> Custom Tracks -> choose ROS_Cfam_1.0 -> paste/upload this file.\n")