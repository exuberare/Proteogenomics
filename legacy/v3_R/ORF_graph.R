# ─────────────────────────────────────────────────────────────────────────────
# ORF-unique peptide evidence across five tissues  (slide "Results – Detection")
# Faceted bar chart: one panel per metric, free x-scales, coloured by prep.
# install.packages(c("ggplot2","dplyr","tidyr","forcats","tibble"))  # if needed
# ─────────────────────────────────────────────────────────────────────────────
library(ggplot2); library(dplyr); library(tidyr); library(forcats); library(tibble)

det <- tribble(
  ~Tissue,  ~Prep,    ~PSMs, ~Peptides, ~CandidateORFs,
  "Brain",  "S-Trap",  243,   33,        21,
  "Heart",  "S-Trap",  483,   36,        15,
  "Liver",  "S-Trap",  112,   12,         7,
  "Lung",   "PAC",     236,   26,        17,
  "Kidney", "PAC",     465,   31,        19
)

# order tissues by ORF-unique PSMs (highest on top), consistent across all panels
lev <- det |> arrange(PSMs) |> pull(Tissue)

long <- det |>
  mutate(Tissue = factor(Tissue, levels = lev)) |>
  pivot_longer(c(PSMs, Peptides, CandidateORFs),
               names_to = "Metric", values_to = "Value") |>
  mutate(Metric = recode_factor(Metric,
                                PSMs          = "ORF-unique PSMs",
                                Peptides      = "Distinct peptides",
                                CandidateORFs = "Candidate ORFs"))

pal <- c("S-Trap" = "#2C6E7F", "PAC" = "#B5651D")   # matches your deck palette

p <- ggplot(long, aes(Value, Tissue, fill = Prep)) +
  geom_col(width = 0.68) +
  geom_text(aes(label = Value), hjust = -0.18, size = 3.6,
            fontface = "bold", colour = "#33383D") +
  facet_wrap(~ Metric, scales = "free_x", nrow = 1) +
  scale_fill_manual(values = pal, name = "Preparation") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(
    title    = "ORF-unique peptide evidence across five tissues",
    subtitle = expression(italic("Canis lupus familiaris")~"· PRIDE PXD052932"),
    caption  = "S-Trap: brain, heart, liver · PAC: lung, kidney — preparation differs, so compare counts within a prep.",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(colour = "grey90"),
    panel.spacing      = unit(1.1, "lines"),
    strip.text    = element_text(face = "bold", size = 13, colour = "#14293A"),
    plot.title    = element_text(face = "bold", size = 16, colour = "#14293A"),
    plot.subtitle = element_text(colour = "grey35", margin = margin(b = 8)),
    plot.caption  = element_text(colour = "grey45", hjust = 0, size = 9),
    legend.position = "top", legend.justification = "left",
    axis.text.y = element_text(face = "bold", colour = "#33383D"),
    axis.text.x = element_text(colour = "grey40")
  )

# white background for a light slide; use bg = "transparent" to drop onto any slide
ggsave("ORF_detection_by_tissue.png", p, width = 11, height = 4.6, dpi = 300, bg = "white")

# ── Alternative: single punchy lollipop of the headline metric (ORF-unique PSMs)
# p2 <- ggplot(det, aes(PSMs, fct_reorder(Tissue, PSMs), colour = Prep)) +
#   geom_segment(aes(x = 0, xend = PSMs, yend = Tissue), linewidth = 1.2) +
#   geom_point(size = 5) +
#   geom_text(aes(label = PSMs), hjust = -0.5, fontface = "bold", size = 4) +
#   scale_colour_manual(values = pal, name = "Preparation") +
#   scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
#   labs(title = "ORF-unique PSMs per tissue", x = "ORF-unique PSMs", y = NULL) +
#   theme_minimal(base_size = 13)
# ggsave("ORF_PSMs_lollipop.png", p2, width = 8, height = 4, dpi = 300, bg = "white")