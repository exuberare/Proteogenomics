# ═══════════════════════════════════════════════════════════════════════════
#  Canine dark-proteome — publication-grade figures  (v3, final)
#  install.packages(c("ggplot2","dplyr","tidyr","readr","readxl","ggridges",
#                     "ggrepel","scales","forcats","pheatmap","UpSetR"))
# ═══════════════════════════════════════════════════════════════════════════
library(ggplot2); library(dplyr); library(tidyr); library(readr)
library(readxl);  library(forcats); library(scales)

PROJECT <- "D:/dark_proteome_final"
META    <- file.path(PROJECT, "orfs/lncRNA_orf_metadata_20aa.csv")
EVID    <- file.path(PROJECT, "reports/dark_proteome_20aa_evidence.xlsx")

OUT <- file.path(PROJECT, "figures"); dir.create(OUT, showWarnings = FALSE); setwd(OUT)

TEAL <- "#2C6E7F"; AMBER <- "#B5651D"; NAVY <- "#14293A"; GREY <- "#8CA0AD"; ink <- "#33383D"
base_theme <- theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        plot.title    = element_text(face = "bold", colour = NAVY, size = 16),
        plot.subtitle = element_text(colour = "grey35", margin = margin(b = 8)),
        plot.caption  = element_text(colour = "grey45", hjust = 0, size = 9),
        axis.title    = element_text(colour = ink))

meta <- read_csv(META, show_col_types = FALSE)
ev   <- read_excel(EVID, sheet = "Evidence_matrix")
psm  <- tryCatch(read_excel(EVID, sheet = "PSM_by_tissue"), error = function(e) NULL)

detected_ids <- ev$orf_id
cand_ids     <- ev$orf_id[grepl("CANDIDATE", ev$final_verdict, ignore.case = TRUE)]

# ───────────────────────────────────────────────────────────────────────────
# 1. RIDGELINE — ORF length: all predicted vs detected vs candidate
# ───────────────────────────────────────────────────────────────────────────
library(ggridges)
len <- bind_rows(
  tibble(aa = meta$aa_length,                                Set = "All predicted ORFs"),
  tibble(aa = meta$aa_length[meta$orf_id %in% detected_ids], Set = "Detected (52)"),
  tibble(aa = meta$aa_length[meta$orf_id %in% cand_ids],     Set = "Candidates (15)")
) |> mutate(Set = factor(Set, levels = c("All predicted ORFs","Detected (52)","Candidates (15)")))

p1 <- ggplot(len, aes(x = aa, y = Set, fill = Set)) +
  geom_density_ridges(scale = 1.15, alpha = 0.85, colour = "white", linewidth = 0.3,
                      quantile_lines = TRUE, quantiles = 2) +
  scale_x_log10(labels = label_number(), breaks = c(20,50,100,300,1000)) +
  scale_fill_manual(values = c(GREY, TEAL, AMBER), guide = "none") +
  labs(title = "Detected ORFs skew longer than the predicted pool",
       subtitle = "ORF length distribution (log scale); dashed line = median",
       x = "ORF length (aa)", y = NULL, caption = "lncRNA-derived ORFs, \u226520 aa.") +
  base_theme
ggsave("fig1_orf_length_ridgeline.png", p1, width = 9, height = 5, dpi = 300, bg = "white")

# ───────────────────────────────────────────────────────────────────────────
# 2. UPSET — cross-tissue detection of the 52 ORFs  (UpSetR = clean classic)
# ───────────────────────────────────────────────────────────────────────────
library(UpSetR)
runcols <- setdiff(grep("_runs$", names(ev), value = TRUE), "max_runs")
bin <- ev |>
  transmute(across(all_of(runcols), ~ as.integer(.x > 0))) |>
  rename_with(~ sub("_runs$", "", .x))
tissues <- names(bin)

png("fig2_tissue_upset.png", width = 2600, height = 1500, res = 300)
print(upset(as.data.frame(bin), sets = tissues, order.by = "freq",
            main.bar.color = TEAL, sets.bar.color = AMBER, matrix.color = NAVY,
            point.size = 3, line.size = 1, text.scale = 1.4,
            mainbar.y.label = "ORFs", sets.x.label = "ORFs per tissue"))
dev.off()

# ───────────────────────────────────────────────────────────────────────────
# 3. HEATMAP — candidate ORFs x tissue (PSM counts, log1p) — pheatmap
# ───────────────────────────────────────────────────────────────────────────
library(pheatmap)
stopifnot(!is.null(psm))
pcols    <- grep("^PSMs_", names(psm), value = TRUE)
show_ids <- if (length(cand_ids)) cand_ids else detected_ids
M <- as.matrix(psm[match(show_ids, psm$orf_id), pcols]); rownames(M) <- show_ids
colnames(M) <- sub("^PSMs_", "", pcols); Ml <- log1p(M)
ann <- data.frame(Verdict = ev$final_verdict[match(show_ids, ev$orf_id)]); rownames(ann) <- show_ids
ann_cols <- list(Verdict = c("CANDIDATE (non-coding both)" = TEAL, "CANDIDATE" = TEAL,
                             "KNOWN" = GREY, "REJECT (coding in either)" = AMBER))
pheatmap(Ml, color = colorRampPalette(c("#F3F6F7", TEAL))(50),
         annotation_row = ann, annotation_colors = ann_cols, border_color = "white",
         fontsize_row = 8, fontsize_col = 11, main = "ORF-unique PSMs per tissue (candidates)",
         filename = "fig3_candidate_heatmap.png", width = 8, height = 7)

# ───────────────────────────────────────────────────────────────────────────
# 4. BUBBLE "candidate landscape" — robust verdict colours, strong labels only
# ───────────────────────────────────────────────────────────────────────────
library(ggrepel)
tot  <- if (!is.null(psm)) psm |> transmute(orf_id, total_PSMs) else
  ev |> transmute(orf_id, total_PSMs = max_runs)
land <- ev |> left_join(tot, by = "orf_id") |>
  mutate(Verdict = case_when(
    grepl("CANDIDATE", final_verdict, ignore.case = TRUE) ~ "Candidate",
    grepl("REJECT",    final_verdict, ignore.case = TRUE) ~ "Reject",
    grepl("KNOWN",     final_verdict, ignore.case = TRUE) ~ "Known",
    TRUE ~ "Other"),
    Verdict = factor(Verdict, levels = c("Candidate","Known","Reject","Other")))
p4 <- ggplot(land, aes(aa_length, total_PSMs, size = n_tissues, fill = Verdict)) +
  geom_point(shape = 21, colour = "white", alpha = 0.9, stroke = 0.5) +
  geom_text_repel(data = subset(land, Verdict == "Candidate" & (n_tissues >= 2 | total_PSMs >= 3)),
                  aes(label = orf_id), size = 3.2, fontface = "bold", max.overlaps = Inf,
                  seed = 1, min.segment.length = 0, box.padding = 0.6, segment.colour = "grey60") +
  scale_y_log10() + scale_size_continuous(range = c(3, 11), breaks = 1:5) +
  scale_fill_manual(values = c(Candidate = TEAL, Known = GREY, Reject = AMBER, Other = "grey75"),
                    breaks = c("Candidate","Known","Reject"), drop = FALSE) +
  guides(fill = guide_legend(override.aes = list(size = 5))) +
  labs(title = "Candidate landscape",
       subtitle = "Each point an ORF; bubble size = tissues detected",
       x = "ORF length (aa)", y = "Total ORF-unique PSMs (log)", size = "Tissues", fill = NULL) +
  base_theme
ggsave("fig4_candidate_landscape.png", p4, width = 9, height = 6, dpi = 300, bg = "white")

# ───────────────────────────────────────────────────────────────────────────
# 5. CONSERVATION scatter — columns: qseqid sseqid pident length EVALUE bitscore stitle
# ───────────────────────────────────────────────────────────────────────────
cons_path <- file.path(PROJECT, "conservation/conservation_hits.tsv")
if (file.exists(cons_path)) {
  rd <- read.table(cons_path, sep = "\t", header = FALSE, fill = TRUE,
                   quote = "", comment.char = "", stringsAsFactors = FALSE)
  cn <- c("orf","subj","pident","length","evalue","bitscore","stitle")
  names(rd)[seq_len(min(ncol(rd), length(cn)))] <- cn[seq_len(min(ncol(rd), length(cn)))]
  
  cons <- tibble(orf    = as.character(rd$orf),
                 pident = suppressWarnings(as.numeric(rd$pident)),
                 evalue = suppressWarnings(as.numeric(rd$evalue))) |>
    filter(!is.na(pident), !is.na(evalue)) |>
    group_by(orf) |> slice_min(evalue, n = 1, with_ties = FALSE) |> ungroup() |>
    mutate(neglogE = -log10(pmax(evalue, 1e-200)),
           Homolog = factor(ifelse(pident >= 70, "Conserved homolog (\u226570%)",
                                   "Distant homolog (<70%)")))
  p5 <- ggplot(cons, aes(pident, neglogE, colour = Homolog)) +
    geom_point(size = 3.2, alpha = 0.9) +
    geom_vline(xintercept = 70, linetype = 2, colour = GREY) +
    scale_colour_manual(values = c("Conserved homolog (\u226570%)" = AMBER,
                                   "Distant homolog (<70%)" = TEAL)) +
    labs(title = "Conservation of detected ORFs vs Swiss-Prot",
         subtitle = "Only ORFs with a Swiss-Prot hit shown; lineage-specific candidates have none",
         x = "% identity to best hit", y = expression(-log[10]~"(E-value)"),
         colour = NULL, caption = "Vertical guide: 70% identity.") +
    base_theme
  ggsave("fig5_conservation_scatter.png", p5, width = 8.5, height = 6, dpi = 300, bg = "white")
} else message("conservation_hits.tsv not found at: ", cons_path)

cat("Done — figures written to", getwd(), "\n")