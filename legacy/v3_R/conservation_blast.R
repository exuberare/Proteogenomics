###############################################################################
##  CROSS-SPECIES CONSERVATION — candidate dark ORFs vs SwissProt (all species)
##  A novel-in-dog microprotein conserved in human/mouse/other is strong
##  evidence it is real and functional. Reports top ortholog, species,
##  %identity, e-value -> feeds the human/mouse Ribo-seq lookup (GWIPS-viz).
###############################################################################
setwd("D:/dark_proteome_final")
suppressPackageStartupMessages({ library(Biostrings); library(openxlsx) })
dir.create("conservation", showWarnings=FALSE); dir.create("conservation/blastdb", showWarnings=FALSE)
options(timeout=3600)

META_CSV <- "orfs/lncRNA_orf_metadata_20aa.csv"
EVIDENCE <- "reports/dark_proteome_20aa_evidence.xlsx"
SWISSPROT_FA <- "conservation/uniprot_sprot.fasta"   # all-species SwissProt (download below)

## ---- which ORFs to test: your candidates (all in evidence matrix) ---------
## (optionally restrict to novel ones; here we test all candidates)
ev <- read.xlsx(EVIDENCE, sheet="Evidence_matrix")
meta <- read.csv(META_CSV, stringsAsFactors=FALSE)
orf_ids  <- unique(ev$orf_id)
orf_seqs <- setNames(meta$aa_sequence[match(orf_ids, meta$orf_id)], orf_ids)
qfa <- AAStringSet(orf_seqs); names(qfa) <- orf_ids
writeXStringSet(qfa, "conservation/candidate_orfs.fasta")
cat("candidates to test for conservation:", length(orf_ids), "\n")

## ---- SwissProt database (all species, reviewed) ---------------------------
## Download ONCE (small, ~90 MB). If missing, this downloads + unzips it.
if (!file.exists(SWISSPROT_FA)) {
  gz <- paste0(SWISSPROT_FA, ".gz")
  if (!file.exists(gz)) {
    cat("downloading SwissProt (all species, reviewed) ...\n")
    download.file("https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz",
                  gz, mode="wb", quiet=TRUE)
  }
  if (!requireNamespace("R.utils", quietly=TRUE)) install.packages("R.utils")
  R.utils::gunzip(gz, remove=FALSE, overwrite=TRUE)
}
stopifnot("SwissProt FASTA missing" = file.exists(SWISSPROT_FA))

## ---- build BLAST db + run blastp ------------------------------------------
blastp <- Sys.which("blastp"); makeblastdb <- Sys.which("makeblastdb")
stopifnot("BLAST+ not on PATH" = nzchar(blastp) && nzchar(makeblastdb))
if (!file.exists("conservation/blastdb/sprot.phr")) {
  cat("building SwissProt BLAST db ...\n")
  system2(makeblastdb, c("-in", shQuote(SWISSPROT_FA), "-dbtype","prot","-out","conservation/blastdb/sprot"))
}
cat("running blastp vs SwissProt (all species) ...\n")
ofmt <- "6 qseqid sseqid pident length evalue bitscore stitle"
system2(blastp, c("-query","conservation/candidate_orfs.fasta","-db","conservation/blastdb/sprot",
                  "-evalue","1e-3","-max_target_seqs","5",
                  "-outfmt", shQuote(ofmt), "-out","conservation/conservation_hits.tsv"))

bl <- tryCatch(read.delim("conservation/conservation_hits.tsv", header=FALSE, quote="", sep="\t"),
               error=function(e) NULL)

## ---- best hit per ORF + species extraction --------------------------------
cons <- data.frame(orf_id=orf_ids, conserved="no", top_hit=NA_character_,
                   species=NA_character_, pident=NA_real_, evalue=NA_real_,
                   stringsAsFactors=FALSE)
if (!is.null(bl) && nrow(bl)) {
  names(bl) <- c("qseqid","sseqid","pident","length","evalue","bitscore","stitle")
  best <- bl[order(bl$qseqid, -bl$bitscore),]; best <- best[!duplicated(best$qseqid),]
  # species is in the title as "OS=Homo sapiens ..."
  sp <- sub(".*OS=([A-Za-z]+ [a-z]+).*", "\\1", best$stitle); sp[sp==best$stitle] <- NA
  m <- match(cons$orf_id, best$qseqid)
  cons$top_hit <- best$stitle[m]
  cons$species <- sp[m]
  cons$pident  <- best$pident[m]
  cons$evalue  <- best$evalue[m]
  cons$conserved <- ifelse(!is.na(cons$pident), "yes", "no")
}
n_cons <- sum(cons$conserved=="yes")
cat("\n=== CONSERVATION ===\n")
cat("candidates with a SwissProt hit (conserved):", n_cons, "of", nrow(cons), "\n")
cat("species of top hits:\n"); print(sort(table(cons$species), decreasing=TRUE))
print(cons[cons$conserved=="yes", c("orf_id","species","pident","evalue","top_hit")], row.names=FALSE)

## ---- merge into evidence + export -----------------------------------------
m <- match(ev$orf_id, cons$orf_id)
ev$conserved        <- cons$conserved[m]
ev$conserved_species<- cons$species[m]
ev$conserved_pident <- cons$pident[m]
ev$conserved_hit    <- cons$top_hit[m]
wb <- loadWorkbook(EVIDENCE)
if ("Evidence_matrix" %in% names(wb)) removeWorksheet(wb, "Evidence_matrix")
addWorksheet(wb, "Evidence_matrix"); writeData(wb, "Evidence_matrix", ev)
saveWorkbook(wb, EVIDENCE, overwrite=TRUE)
write.xlsx(cons, "conservation/conservation_summary.xlsx", overwrite=TRUE)
cat("\nUpdated evidence matrix + wrote conservation/conservation_summary.xlsx\n")
cat("NEXT: for any hit in Homo sapiens / Mus musculus, look up that ortholog in\n")
cat("      GWIPS-viz (human/mouse) to check for ribosome-footprint (translation) signal.\n")