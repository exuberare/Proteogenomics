###############################################################################
##  CANINE DARK PROTEOME — GENERAL ANALYSIS PIPELINE
##  Tissue-agnostic: list your tissues + their FragPipe output folders ONCE
##  in the CONFIG block below. Add a new tissue later = add one line. Nothing
##  else in the script changes.
##
##  Runs, in order:
##    1. Pool PSMs per tissue -> ORF-unique peptides -> per-tissue candidates
##    2. Cross-tissue table (which ORF in which tissues; reproducibility)
##    3. BLAST novelty (peptide exact-match + protein BLASTp)
##    4. tBLASTn genomic mapping (+ CDS-overlap flag)
##    5. One consolidated evidence matrix exported to Excel
###############################################################################

## ======================= CONFIG — edit only this =========================
PROJECT  <- "D:/dark_proteome_final"          # working project root
META_CSV <- "orfs/lncRNA_orf_metadata_20aa.csv"
GENOME_FA<- "reference_genome/genome/Canis_lupus_familiaris.ROS_Cfam_1.0.dna.toplevel.fa"
GTF_PATH <- "reference_genome/genome/Canis_lupus_familiaris.ROS_Cfam_1.0.115.gtf"
PEP_FA   <- "reference_genome/genome/Canis_lupus_familiaris.ROS_Cfam_1.0.pep.all.fa"

# ONE line per tissue: name = FragPipe output folder. Add more any time.
# Commas after every line EXCEPT the last. No '#' in front of active tissues.
TISSUES <- c(
  brain  = "D:/canine_lncRNA_proteogenomics/proteomics/fragpipe_brain_out_20aa",
  heart  = "D:/canine_lncRNA_proteogenomics/proteomics/fragpipe_heart_out_20aa",
  liver  = "D:/canine_lncRNA_proteogenomics/proteomics/fragpipe_liver_out_20aa",
  lung   = "D:/canine_lncRNA_proteogenomics/proteomics/fragpipe_lung_out_20aa",
  kidney = "D:/canine_lncRNA_proteogenomics/proteomics/fragpipe_kidney_out_20aa"
)

MIN_RUNS_REPRO <- 2          # runs within a tissue to call "reproducible"
RUN_NOVELTY    <- TRUE       # needs UniProt+RefSeq in novelty/ (auto-copied)
RUN_TBLASTN    <- TRUE       # needs BLAST+ on PATH and the genome FASTA
OUT_TAG        <- "20aa"     # suffix for output files
## =========================================================================

setwd(PROJECT)
suppressPackageStartupMessages({
  library(Biostrings); library(openxlsx)
  if (RUN_TBLASTN) { library(rtracklayer); library(GenomicRanges) }
})
for (d in c("reports","novelty","novelty/blast","tblastn","tblastn/blastdb"))
  dir.create(d, recursive=TRUE, showWarnings=FALSE)
meta <- read.csv(META_CSV, stringsAsFactors=FALSE)

## ---- helper: pool a tissue's PSMs -> ORF-unique --------------------------
read_tissue <- function(folder){
  pf <- list.files(folder, pattern="psm\\.tsv$", recursive=TRUE, full.names=TRUE)
  if (length(pf)==0) { warning("no psm.tsv in ", folder); return(NULL) }
  psm <- do.call(rbind, lapply(pf, function(f){
    d <- read.delim(f, check.names=FALSE); d$run <- basename(dirname(f)); d }))
  psm$orf_lead <- grepl("DARK_ORF", psm$"Protein")
  psm$has_ref  <- grepl("ENSCAFP", psm$"Protein") | grepl("ENSCAFP", psm$"Mapped Proteins")
  psm[psm$orf_lead & !psm$has_ref, , drop=FALSE]
}

## ---- 1. per-tissue candidates --------------------------------------------
cat("=== Step 1: reading tissues ===\n")
uniq_list <- list(); all_uniq <- NULL
for (t in names(TISSUES)) {
  u <- read_tissue(TISSUES[[t]])
  if (is.null(u)) { cat("  ", t, ": FOLDER NOT READY — skipped\n"); next }
  u$tissue <- t; uniq_list[[t]] <- u; all_uniq <- rbind(all_uniq, u)
  cat(sprintf("  %-7s: %d ORF-unique PSMs | %d peptides | %d ORFs\n",
              t, nrow(u), length(unique(u$Peptide)), length(unique(u$"Protein"))))
}
tissues_done <- names(uniq_list)
stopifnot("No tissue folders ready." = length(tissues_done) > 0)

## per-tissue per-ORF summary
summarise_tissue <- function(u){
  ids <- unique(u$"Protein")
  data.frame(orf_id=ids,
             n_pep=as.integer(tapply(u$Peptide,u$"Protein",function(x)length(unique(x)))[ids]),
             n_run=as.integer(tapply(u$run,u$"Protein",function(x)length(unique(x)))[ids]),
             n_psm=as.integer(table(u$"Protein")[ids]),
             peptides=tapply(u$Peptide,u$"Protein",function(x)paste(sort(unique(x)),collapse="; "))[ids],
             stringsAsFactors=FALSE)
}
per_tissue <- lapply(uniq_list, summarise_tissue)

## ---- 2. cross-tissue matrix ----------------------------------------------
cat("\n=== Step 2: cross-tissue ===\n")
all_ids <- sort(unique(all_uniq$"Protein"))
X <- data.frame(orf_id=all_ids,
                aa_length = meta$aa_length[match(all_ids, meta$orf_id)],
                num_exons = meta$num_exons[match(all_ids, meta$orf_id)],
                transcript_id = meta$transcript_id[match(all_ids, meta$orf_id)],
                gene_id = meta$gene_id[match(all_ids, meta$orf_id)],
                stringsAsFactors=FALSE)
for (t in tissues_done) {
  st <- per_tissue[[t]]
  X[[paste0(t,"_runs")]] <- st$n_run[match(all_ids, st$orf_id)]
  X[[paste0(t,"_pep")]]  <- st$n_pep[match(all_ids, st$orf_id)]
}
runcols <- paste0(tissues_done,"_runs")
X[runcols][is.na(X[runcols])] <- 0
X[paste0(tissues_done,"_pep")][is.na(X[paste0(tissues_done,"_pep")])] <- 0
X$n_tissues <- rowSums(X[runcols] > 0)
X$max_runs  <- apply(X[runcols], 1, max)
X$reproducible <- X$max_runs >= MIN_RUNS_REPRO   # in >=1 tissue
X <- X[order(-X$n_tissues, -X$max_runs, -X$aa_length), ]
cat("  total candidates:", nrow(X), "| in >=2 tissues:", sum(X$n_tissues>=2),
    "| reproducible in >=1 tissue:", sum(X$reproducible), "\n")

## ---- 3. BLAST novelty -----------------------------------------------------
X$peptide_novel <- NA; X$blast_pident <- NA_real_; X$blast_evalue <- NA_real_; X$verdict <- NA
if (RUN_NOVELTY) {
  cat("\n=== Step 3: novelty ===\n")
  cpf <- function(f, src) if(!file.exists(f)&&file.exists(src)) file.copy(src,f)
  cpf("novelty/uniprot_dog_9615.fasta","D:/canine_r_msc_project/novelty/uniprot_dog_9615.fasta")
  cpf("novelty/refseq_dog_protein.faa","D:/canine_r_msc_project/novelty/refseq_dog_protein.faa")
  kf <- c("novelty/uniprot_dog_9615.fasta","novelty/refseq_dog_protein.faa", PEP_FA)
  kf <- kf[file.exists(kf)]
  if (length(kf) >= 2) {
    known <- as(Reduce(append, lapply(kf, readAAStringSet)), "AAStringSet")
    il <- function(x) chartr("I","L",toupper(x)); ref_blob <- paste(il(as.character(known)),collapse="*")
    # peptide-level, I/L-aware, pooled across tissues
    pep_tbl <- unique(all_uniq[,c("Protein","Peptide")])
    peps <- sort(unique(pep_tbl$Peptide))
    found <- vapply(il(peps), function(p) grepl(p, ref_blob, fixed=TRUE), logical(1))
    pep_tbl$found <- found[match(pep_tbl$Peptide, peps)]
    pn <- !tapply(pep_tbl$found, pep_tbl$Protein, any)
    X$peptide_novel <- pn[X$orf_id]
    # protein BLASTp
    blastp <- Sys.which("blastp"); makeblastdb <- Sys.which("makeblastdb")
    if (nzchar(blastp) && nzchar(makeblastdb)) {
      qfa <- AAStringSet(setNames(meta$aa_sequence[match(X$orf_id,meta$orf_id)], X$orf_id))
      writeXStringSet(qfa, "novelty/candidate_orfs.fasta")
      if (!file.exists("novelty/blast/known_dog.phr")) {
        writeXStringSet(known,"novelty/known_dog_combined.fasta")
        system2(makeblastdb,c("-in","novelty/known_dog_combined.fasta","-dbtype","prot","-out","novelty/blast/known_dog")) }
      system2(blastp, c("-query","novelty/candidate_orfs.fasta","-db","novelty/blast/known_dog",
                        "-evalue","1e-2","-max_target_seqs","5","-num_threads","4",
                        "-outfmt",shQuote("6 qseqid sseqid pident length qlen slen evalue bitscore"),
                        "-out","novelty/blast/orf_hits.tsv"))
      bl <- tryCatch(read.delim("novelty/blast/orf_hits.tsv",header=FALSE),error=function(e)NULL)
      if (!is.null(bl) && nrow(bl)) {
        names(bl)<-c("q","s","pid","len","qlen","slen","ev","bit"); bl$cov<-100*bl$len/bl$qlen
        bb<-bl[order(bl$q,-bl$bit),]; bb<-bb[!duplicated(bb$q),]; m<-match(X$orf_id,bb$q)
        X$blast_pident<-bb$pid[m]; X$blast_evalue<-bb$ev[m]; X$blast_cov<-round(bb$cov[m],1)
      }
    } else cat("  BLAST+ not found — peptide-level novelty only.\n")
    vf <- function(i){ p<-X$peptide_novel[i]; pid<-X$blast_pident[i]; ev<-X$blast_evalue[i]
    cov<-if(!is.null(X$blast_cov)) X$blast_cov[i] else NA
    kf2<- !is.na(pid)&&pid>=98&&!is.na(cov)&&cov>=90; ho<- !is.na(ev)&&ev<=1e-5&&!kf2
    if(is.na(p)) NA else if(!p||kf2) "KNOWN" else if(ho) "NOVEL (conserved homolog)"
    else if(is.na(pid)) "NOVEL (no BLAST hit)" else "NOVEL (no significant hit)" }
    X$verdict <- vapply(seq_len(nrow(X)), vf, character(1))
    cat("  verdicts:\n"); print(table(X$verdict, useNA="ifany"))
  } else cat("  UniProt/RefSeq not available — novelty skipped.\n")
}

## ---- 4. tBLASTn genomic mapping ------------------------------------------
X$chromosome<-NA; X$g_start<-NA; X$g_end<-NA; X$strand<-NA; X$overlaps_known_CDS<-NA
if (RUN_TBLASTN) {
  cat("\n=== Step 4: tBLASTn ===\n")
  tblastn<-Sys.which("tblastn"); makeblastdb<-Sys.which("makeblastdb")
  if (nzchar(tblastn) && nzchar(makeblastdb) && file.exists(GENOME_FA)) {
    qfa <- AAStringSet(setNames(meta$aa_sequence[match(X$orf_id,meta$orf_id)], X$orf_id))
    writeXStringSet(qfa, "tblastn/candidate_orfs.fasta")
    if (!file.exists("tblastn/blastdb/genome.nsq") && !file.exists("tblastn/blastdb/genome.nin")) {
      cat("  building genome DB (one-time) ...\n")
      system2(makeblastdb, c("-in",shQuote(GENOME_FA),"-dbtype","nucl","-out","tblastn/blastdb/genome")) }
    system2(tblastn, c("-query","tblastn/candidate_orfs.fasta","-db","tblastn/blastdb/genome",
                       "-evalue","1e-5","-max_target_seqs","5","-db_gencode","1","-num_threads","4",
                       "-outfmt",shQuote("6 qseqid sseqid pident length qstart qend sstart send evalue bitscore"),
                       "-out","tblastn/orf_genome_hits.tsv"))
    g <- tryCatch(read.delim("tblastn/orf_genome_hits.tsv",header=FALSE),error=function(e)NULL)
    if (!is.null(g) && nrow(g)) {
      names(g)<-c("orf","chr","pid","len","qs","qe","ss","se","ev","bit")
      gb<-g[order(g$orf,-g$bit),]; gb<-gb[!duplicated(gb$orf),]
      gb$strand<-ifelse(gb$ss<=gb$se,"+","-"); gb$gs<-pmin(gb$ss,gb$se); gb$ge<-pmax(gb$ss,gb$se)
      m<-match(X$orf_id,gb$orf)
      X$chromosome<-gb$chr[m]; X$g_start<-gb$gs[m]; X$g_end<-gb$ge[m]; X$strand<-gb$strand[m]
      gtf<-rtracklayer::import(GTF_PATH); cds<-gtf[gtf$type=="CDS"]
      ok<-!is.na(X$chromosome)
      if (any(ok)) {
        hit<-GRanges(X$chromosome[ok],IRanges(X$g_start[ok],X$g_end[ok]),strand=X$strand[ok])
        ov<-findOverlaps(hit,cds,ignore.strand=TRUE)
        flag<-rep(FALSE,sum(ok)); flag[queryHits(ov)]<-TRUE; X$overlaps_known_CDS[ok]<-flag
      }
      cat("  mapped:", sum(ok), "of", nrow(X),
          "| overlapping known CDS:", sum(X$overlaps_known_CDS,na.rm=TRUE), "\n")
    }
  } else cat("  BLAST+ or genome FASTA missing — tBLASTn skipped.\n")
}

## ---- 5. consolidated export ----------------------------------------------
cat("\n=== Step 5: export ===\n")
X$locus <- ifelse(is.na(X$chromosome), NA,
                  sprintf("%s:%d-%d(%s)", X$chromosome, X$g_start, X$g_end, X$strand))
print(X[, intersect(c("orf_id","aa_length","n_tissues","max_runs","peptide_novel","verdict","overlaps_known_CDS"), names(X))], row.names=FALSE)

write.xlsx(list(Evidence_matrix = X,
                All_ORF_unique_PSMs = all_uniq),
           sprintf("reports/dark_proteome_%s_evidence.xlsx", OUT_TAG), overwrite=TRUE)
cat("\nWRITTEN reports/dark_proteome_", OUT_TAG, "_evidence.xlsx\n", sep="")
cat("Tissues included:", paste(tissues_done, collapse=", "), "\n")