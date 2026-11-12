###############################################################################
##  NOVELTY CHECK — 20 aa brain candidates  (D:/dark_proteome_final)
##  Stage 1: peptide exact-match vs known dog proteins (I/L-aware)
##  Stage 2: protein-level BLASTp vs UniProt + RefSeq + Ensembl
###############################################################################
setwd("D:/dark_proteome_final")
suppressPackageStartupMessages({ library(Biostrings); library(openxlsx) })
dir.create("novelty", showWarnings=FALSE); dir.create("novelty/blast", showWarnings=FALSE)
options(timeout=3600)

## ---- reference databases (copy from old project if not already here) ------
ensure <- function(f, src){ if(!file.exists(f) && file.exists(src)){ cat("copying",basename(f),"\n"); file.copy(src,f) } }
dir.create("novelty", showWarnings=FALSE)
UNIPROT <- "novelty/uniprot_dog_9615.fasta"
REFSEQ  <- "novelty/refseq_dog_protein.faa"
ENSEMBL <- "reference_genome/genome/Canis_lupus_familiaris.ROS_Cfam_1.0.pep.all.fa"
ensure(UNIPROT, "D:/canine_r_msc_project/novelty/uniprot_dog_9615.fasta")
ensure(REFSEQ,  "D:/canine_r_msc_project/novelty/refseq_dog_protein.faa")
known_files <- c(UniProt=UNIPROT, RefSeq=REFSEQ, Ensembl=ENSEMBL)
have <- setNames(file.exists(known_files), names(known_files)); print(have)
stopifnot(all(have[c("UniProt","RefSeq")]))
known <- Reduce(append, lapply(known_files[have], readAAStringSet)); known <- as(known,"AAStringSet")
cat("known dog proteins:", length(known), "\n")

## ---- candidates (brain 20aa) ----------------------------------------------
CAND <- "reports/brain_20aa_candidates.xlsx"
cand  <- read.xlsx(CAND, sheet="Brain_20aa_candidates")
uniqP <- read.xlsx(CAND, sheet="All_ORF_unique_PSMs")
meta  <- read.csv("orfs/lncRNA_orf_metadata_20aa.csv", stringsAsFactors=FALSE)
orf_ids  <- unique(cand$orf_id)
orf_seqs <- setNames(meta$aa_sequence[match(orf_ids, meta$orf_id)], orf_ids)
pep_tbl  <- unique(uniqP[, c("Protein","Peptide")])
cat("candidates:", length(orf_ids), "| unique peptides:", length(unique(pep_tbl$Peptide)), "\n")

## ---- STAGE 1: peptide exact-match (I/L-aware) -----------------------------
il <- function(x) chartr("I","L",toupper(x))
ref_blob <- paste(il(as.character(known)), collapse="*")
peptides <- sort(unique(pep_tbl$Peptide))
pep_found <- vapply(il(peptides), function(p) grepl(p, ref_blob, fixed=TRUE), logical(1))
pep_res <- data.frame(Peptide=peptides, found_in_known_dog=pep_found, row.names=NULL)
pep_tbl$found <- pep_res$found_in_known_dog[match(pep_tbl$Peptide, pep_res$Peptide)]
orf_peptide_novel <- !tapply(pep_tbl$found, pep_tbl$Protein, any)[orf_ids]
cat("\nSTAGE 1: peptides in a known dog protein:", sum(pep_found), "of", length(pep_found), "\n")
cat("  ORFs peptide-level NOVEL:", sum(orf_peptide_novel), "of", length(orf_ids), "\n")

## ---- STAGE 2: protein BLASTp ----------------------------------------------
blastp <- Sys.which("blastp"); makeblastdb <- Sys.which("makeblastdb")
orf_blast <- data.frame(orf_id=orf_ids, best_pident=NA_real_, best_cov=NA_real_,
                        best_evalue=NA_real_, best_hit=NA_character_, stringsAsFactors=FALSE)
if (nzchar(blastp) && nzchar(makeblastdb)) {
  cat("\nSTAGE 2: running blastp ...\n")
  qfa <- AAStringSet(orf_seqs); names(qfa) <- orf_ids
  writeXStringSet(qfa, "novelty/candidate_orfs.fasta")
  if (!file.exists("novelty/blast/known_dog.phr"))
  { writeXStringSet(known,"novelty/known_dog_combined.fasta")
    system2(makeblastdb, c("-in","novelty/known_dog_combined.fasta","-dbtype","prot","-out","novelty/blast/known_dog")) }
  ofmt <- "6 qseqid sseqid pident length qlen slen evalue bitscore"
  system2(blastp, c("-query","novelty/candidate_orfs.fasta","-db","novelty/blast/known_dog",
                    "-evalue","1e-2","-max_target_seqs","5","-outfmt",shQuote(ofmt),"-out","novelty/blast/orf_hits.tsv"))
  bl <- tryCatch(read.delim("novelty/blast/orf_hits.tsv",header=FALSE), error=function(e) NULL)
  if (!is.null(bl) && nrow(bl)) {
    names(bl) <- c("qseqid","sseqid","pident","length","qlen","slen","evalue","bitscore")
    bl$cov <- 100*bl$length/bl$qlen
    best <- bl[order(bl$qseqid,-bl$bitscore),]; best <- best[!duplicated(best$qseqid),]
    m <- match(orf_ids, best$qseqid)
    orf_blast$best_pident<-best$pident[m]; orf_blast$best_cov<-round(best$cov[m],1)
    orf_blast$best_evalue<-best$evalue[m]; orf_blast$best_hit<-best$sseqid[m]
  }
} else cat("\nSTAGE 2: BLAST+ not on PATH — skipped (Stage 1 still valid).\n")

## ---- verdict + export -----------------------------------------------------
verdict <- function(i){ pn<-orf_peptide_novel[[i]]; pid<-orf_blast$best_pident[i]
cov<-orf_blast$best_cov[i]; ev<-orf_blast$best_evalue[i]
kf<- !is.na(pid)&&pid>=98&&!is.na(cov)&&cov>=90; ho<- !is.na(ev)&&ev<=1e-5&&!kf
if(!pn||kf) "KNOWN (matches known dog protein)"
else if(ho) "NOVEL in dog, conserved homolog"
else if(is.na(pid)) "NOVEL (no BLAST hit)"
else "NOVEL (no significant hit)" }
mm <- function(c) meta[[c]][match(orf_ids, meta$orf_id)]
res <- data.frame(orf_id=orf_ids, aa_length=mm("aa_length"),
                  n_peptides=cand$n_peptides[match(orf_ids,cand$orf_id)],
                  n_runs=cand$n_runs[match(orf_ids,cand$orf_id)],
                  transcript_id=mm("transcript_id"), gene_id=mm("gene_id"),
                  peptide_novel=orf_peptide_novel[orf_ids],
                  blast_pident=orf_blast$best_pident, blast_evalue=orf_blast$best_evalue,
                  stringsAsFactors=FALSE)
res$verdict <- vapply(seq_len(nrow(res)), verdict, character(1))
res <- res[order(res$verdict, -res$n_runs, -res$n_peptides),]
cat("\n=== NOVELTY VERDICTS (brain 20aa) ===\n"); print(table(res$verdict))
print(res[,c("orf_id","aa_length","n_peptides","n_runs","peptide_novel","verdict")], row.names=FALSE)
write.xlsx(list(Novelty_verdict=res, Peptide_matches=pep_res),
           "novelty/brain_20aa_novelty_verdict.xlsx", overwrite=TRUE)
cat("\nwritten novelty/brain_20aa_novelty_verdict.xlsx\n")