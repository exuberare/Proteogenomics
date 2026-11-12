###############################################################################
##  CANINE DARK PROTEOME — 20 aa RE-RUN  (Objective 1, database build)
##  Fresh output in D:/dark_proteome_final ; ORF floor lowered 50 -> 20 aa.
##
##  Reuses the already-extracted lncRNA transcripts from the 50-aa project if
##  present (fast); otherwise rebuilds them from the Ensembl genome.
##  ORF rule unchanged except length: complete ORF (ATG -> stop), >= 20 aa,
##  three forward frames, longest ORF per stop window.
###############################################################################

## ---- setup ----------------------------------------------------------------
setwd("D:/dark_proteome_final")
suppressPackageStartupMessages({
  library(rtracklayer); library(GenomicFeatures); library(txdbmaker)
  library(Biostrings);  library(GenomicRanges);   library(GenomeInfoDb)
  library(openxlsx)
})
for (d in c("reference_genome/genome","transcripts","orfs","databases","reports"))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

MIN_AA  <- 20L                       # <-- lowered from 50
ORF_TAG <- "DARK_ORF"
OLD     <- "D:/canine_r_msc_project" # source of reference files / transcripts

## ---- reference paths (copy from old project if missing) -------------------
GENOME_FA <- "reference_genome/genome/Canis_lupus_familiaris.ROS_Cfam_1.0.dna.toplevel.fa"
GTF_PATH  <- "reference_genome/genome/Canis_lupus_familiaris.ROS_Cfam_1.0.115.gtf"
PEP_FA    <- "reference_genome/genome/Canis_lupus_familiaris.ROS_Cfam_1.0.pep.all.fa"
for (f in c(GENOME_FA, GTF_PATH, PEP_FA)) {
  if (!file.exists(f)) {
    src <- file.path(OLD, f)
    if (file.exists(src)) { cat("copying", basename(f), "...\n"); file.copy(src, f) }
  }
}
stopifnot("GTF missing"     = file.exists(GTF_PATH),
          "proteome missing"= file.exists(PEP_FA))

## ---- lncRNA selection (from GTF) ------------------------------------------
gtf <- rtracklayer::import(GTF_PATH)
tx_rows    <- gtf[gtf$type == "transcript"]
lnc_tx_ids <- unique(tx_rows$transcript_id[tx_rows$gene_biotype == "lncRNA"])
cat("lncRNA transcripts:", length(lnc_tx_ids), "\n")
tx2gene <- setNames(tx_rows$gene_id[match(lnc_tx_ids, tx_rows$transcript_id)], lnc_tx_ids)
gname   <- tx_rows$gene_name[match(lnc_tx_ids, tx_rows$transcript_id)]
gname[is.na(gname) | gname == ""] <- tx2gene[is.na(gname) | gname == ""]
tx2name <- setNames(gname, lnc_tx_ids)

## ---- transcripts: reuse if available, else extract from genome ------------
OLD_TX <- file.path(OLD, "transcripts/lncRNA_transcripts.fa")
if (file.exists("transcripts/lncRNA_transcripts.fa")) {
  tx_seqs <- readDNAStringSet("transcripts/lncRNA_transcripts.fa")
} else if (file.exists(OLD_TX)) {
  cat("reusing transcripts from old project ...\n")
  file.copy(OLD_TX, "transcripts/lncRNA_transcripts.fa")
  tx_seqs <- readDNAStringSet("transcripts/lncRNA_transcripts.fa")
} else {
  cat("extracting transcripts from genome (slow) ...\n")
  stopifnot("genome FASTA missing" = file.exists(GENOME_FA))
  mk_txdb <- if (requireNamespace("txdbmaker", quietly=TRUE))
    txdbmaker::makeTxDbFromGFF else GenomicFeatures::makeTxDbFromGFF
  txdb <- mk_txdb(GTF_PATH, format = "gtf")
  exby <- exonsBy(txdb, by = "tx", use.names = TRUE)
  exby <- exby[names(exby) %in% lnc_tx_ids]
  genome <- readDNAStringSet(GENOME_FA); names(genome) <- sub("\\s.*","",names(genome))
  exby <- keepSeqlevels(exby, intersect(seqlevels(exby), names(genome)), pruning.mode="coarse")
  ex <- unlist(exby, use.names=FALSE); txid <- rep(names(exby), elementNROWS(exby))
  s <- subseq(genome[as.character(seqnames(ex))], start=start(ex), end=end(ex))
  o <- order(match(txid,names(exby)), start(ex)); s<-s[o]; txid<-txid[o]; strd<-as.character(strand(ex))[o]
  spl <- split(s, factor(txid, levels=names(exby)))
  tx_seqs <- DNAStringSet(lapply(spl, function(x) unlist(x)))
  neg <- tapply(strd, factor(txid,levels=names(exby)), function(z) z[1]=="-")[names(tx_seqs)]
  tx_seqs[neg] <- reverseComplement(tx_seqs[neg])
  writeXStringSet(tx_seqs, "transcripts/lncRNA_transcripts.fa")
}
names(tx_seqs) <- sub("\\s.*","",names(tx_seqs))
tx_seqs <- tx_seqs[names(tx_seqs) %in% lnc_tx_ids]
cat("transcripts loaded:", length(tx_seqs), "\n")

## ---- exon counts (for metadata) -------------------------------------------
ex_tx <- gtf[gtf$type=="exon" & gtf$transcript_id %in% names(tx_seqs)]
num_exons_v <- table(factor(ex_tx$transcript_id, levels=names(tx_seqs)))

## ---- ORF caller: complete ORF, ATG->stop, >= MIN_AA, 3 fwd frames ---------
.bases <- c("T","C","A","G")
.aa <- strsplit("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG","")[[1]]
.cod<-character(0); for(a in .bases)for(b in .bases)for(c in .bases).cod<-c(.cod,paste0(a,b,c))
GC <- setNames(.aa,.cod)
translate_frame <- function(nt){ nt<-toupper(gsub("U","T",nt)); k<-nchar(nt)%/%3L
if(k==0L) return(""); idx<-seq(1L,by=3L,length.out=k)
out<-GC[substring(nt,idx,idx+2L)]; out[is.na(out)]<-"X"; paste(out,collapse="") }
find_orfs <- function(seq,min_aa=MIN_AA){ L<-nchar(seq); rows<-list()
for(frame in 0:2){ prot<-strsplit(translate_frame(substring(seq,frame+1L,L)),"")[[1]]
np<-length(prot); start<-1L; i<-1L
while(i<=np){ if(prot[i]=="*"){ win<-if(i>start)prot[start:(i-1L)]else character(0)
if(length(win)&&any(win=="M")){ m<-which(win=="M")[1L]; aalen<-length(win)-m+1L
if(aalen>=min_aa){ p0<-start+m-1L; nt_start<-frame+(p0-1L)*3L+1L
rows[[length(rows)+1L]]<-data.frame(frame=frame+1L,nt_start=nt_start,
                                    nt_end=nt_start+3L*aalen+2L,aa_length=aalen,
                                    aa_sequence=paste(win[m:length(win)],collapse=""),stringsAsFactors=FALSE) } }
start<-i+1L }; i<-i+1L } }
if(!length(rows))return(NULL); df<-do.call(rbind,rows); df[order(df$frame,df$nt_start),] }

## ---- run ORF prediction ---------------------------------------------------
seqs_chr <- as.character(tx_seqs); parts <- vector("list",length(seqs_chr))
for(i in seq_along(seqs_chr)){ o<-find_orfs(seqs_chr[i])
if(!is.null(o)){ tid<-names(seqs_chr)[i]; o$transcript_id<-tid
o$gene_id<-unname(tx2gene[tid]); o$gene_name<-unname(tx2name[tid])
o$num_exons<-as.integer(num_exons_v[tid]); parts[[i]]<-o }
if(i%%1000==0) cat("  ORF scan:",i,"transcripts\n") }
meta <- do.call(rbind,parts)
meta <- meta[order(match(meta$transcript_id,names(tx_seqs)),meta$frame,meta$nt_start),]
meta$orf_id <- sprintf("%s_%06d",ORF_TAG,seq_len(nrow(meta)))
meta <- meta[,c("orf_id","transcript_id","gene_id","gene_name","frame",
                "nt_start","nt_end","aa_length","num_exons","aa_sequence")]
cat("\ncandidate ORFs (>=",MIN_AA,"aa):", nrow(meta), "\n")
cat("  length aa  median/mean/max:", median(meta$aa_length),"/",
    round(mean(meta$aa_length),1),"/",max(meta$aa_length),"\n")
cat("  ORFs in 20-49 aa range (new):", sum(meta$aa_length<50), "\n")

# integrity check
stopifnot(sum(substr(meta$aa_sequence,1,1)!="M")==0,
          sum(grepl("\\*|X",meta$aa_sequence))==0,
          sum(meta$aa_length<MIN_AA)==0)

## ---- write ORF FASTA + metadata -------------------------------------------
hdr <- sprintf("%s transcript=%s gene=%s frame=%d start=%d end=%d length=%daa",
               meta$orf_id,meta$transcript_id,meta$gene_id,meta$frame,meta$nt_start,meta$nt_end,meta$aa_length)
orf_aa <- AAStringSet(meta$aa_sequence); names(orf_aa)<-hdr
writeXStringSet(orf_aa,"orfs/lncRNA_ORFs_20aa.fasta")
write.csv(meta,"orfs/lncRNA_orf_metadata_20aa.csv",row.names=FALSE)

## ---- build combined target database ---------------------------------------
ref <- readAAStringSet(PEP_FA); names(ref)<-sub("\\s.*","",names(ref))
orf_dedup <- orf_aa[!duplicated(as.character(orf_aa))]; names(orf_dedup)<-sub("\\s.*","",names(orf_dedup))
target <- c(ref,orf_dedup); target <- target[!duplicated(names(target))]
writeXStringSet(target,"databases/canine_target_db_20aa.fasta")
cat("\nreference proteins:",length(ref),"| ORFs after dedup:",length(orf_dedup),
    "| TARGET DB:",length(target),"\n")
cat("  -> databases/canine_target_db_20aa.fasta  (load this in FragPipe)\n")

## ---- summary --------------------------------------------------------------
summ <- data.frame(metric=c("lncRNA transcripts","candidate ORFs (>=20 aa)",
                            "new ORFs 20-49 aa","ORFs after dedup","reference proteins","target DB total",
                            "median aa","max aa"),
                   value=c(length(tx_seqs),nrow(meta),sum(meta$aa_length<50),length(orf_dedup),
                           length(ref),length(target),median(meta$aa_length),max(meta$aa_length)))
write.xlsx(list(Summary=summ, All_ORFs=meta),
           "reports/lncRNA_ORFs_20aa_review.xlsx", overwrite=TRUE)
saveRDS(meta,"reports/orf_metadata_20aa.rds")
print(summ)
cat("\nDONE. FragPipe DB: databases/canine_target_db_20aa.fasta\n")