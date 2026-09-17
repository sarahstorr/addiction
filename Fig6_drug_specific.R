
# ---- 1. Settings ----------------------------------------------------------

analysis_dir <- "C:/Users/mrzsjs1/Desktop/Addiction/analysis"
out_dir      <- file.path(analysis_dir, "figures", "Fig6")
PHOS_DIR     <- analysis_dir                          # *_phosphoproteomics_significant.xlsx
METAB_XLSX   <- "C:/Users/mrzsjs1/Desktop/Addiction/analysis/metabolomics final.xlsx"

LINES <- c("KNS-42" = "KNS", "SVG-A" = "SVG")           # label -> filename prefix
DRUGS <- c(ethanol = "ethanol", THC = "thc", cocaine = "cocaine",
           morphine = "morphine", psilocin = "psilocin")
SVG_DRUGS <- c("ethanol", "THC", "cocaine", "psilocin")  # SVG-A morphine excluded (flowcell)

# transcript files; KNS-42 ethanol and psilocin use the matched-control reruns
tx_file <- function(prefix, drug) {
  if (prefix == "KNS" && drug %in% c("ethanol", "psilocin"))
    return(sprintf("KNS_%s_DESeq2_transcript_results_CORRECTED.csv", drug))
  sprintf("%s_%s_DESeq2_transcript_results.csv", prefix, drug)
}
phos_file   <- function(drug) sprintf("%s_phosphoproteomics_significant.xlsx", DRUGS[[drug]])
METAB_SHEET <- "Filtered"
## drug group vs its control (THC vs ethanol vehicle); the sheet reports
## control / drug ratios, so log2FC is negated to give drug vs control
METAB_CONTRAST <- list(ethanol = c("E","UC"), THC = c("T","EC"), cocaine = c("C","UC"),
                       morphine = c("M","UC"), psilocin = c("P","UC"))

# significance rule
PADJ <- 0.05; LFC <- 1
NEAR_P <- 0.05; NEAR_LFC <- 0.5             # near-threshold definition
N_PERM_LAYER <- 100000; N_PERM_LINE <- 20000
SEED <- 20260915
SAVE_SUPPLEMENTARY <- TRUE                  # also save the supplementary panels

# palette
COL_UP <- "#E24B4A"; COL_DOWN <- "#378ADD"; COL_NS <- "#B4B2A9"
COL_CONC <- "#1D9E75"; COL_DISC <- "#E24B4A"
LINE_COLS <- c("KNS-42" = "#1D9E75", "SVG-A" = "#7F77DD")
DRUG_COLS <- c(ethanol = "#E69F00", THC = "#1D9E75", cocaine = "#378ADD",
               morphine = "#E24B4A", psilocin = "#7F77DD")


# ---- 2. Packages ----------------------------------------------------------

pkgs <- c("readr", "readxl", "dplyr", "tidyr", "ggplot2", "scales", "patchwork", "ggrepel")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing)
suppressPackageStartupMessages({
  library(readr); library(readxl); library(dplyr); library(tidyr)
  library(ggplot2); library(scales); library(patchwork); library(ggrepel)
})
select <- dplyr::select; filter <- dplyr::filter
rename <- dplyr::rename; mutate <- dplyr::mutate

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
if (!file.exists(METAB_XLSX)) stop("Metabolomics file not found: ", METAB_XLSX)
set.seed(SEED)
thr_lab <- sprintf("padj < %.2g & |log2FC| > %g", PADJ, LFC)

# ---- Text sizes and theme --------------------------------------------------
# Figures are drawn at their printed size, so SZ values are the point sizes
# seen in the journal at 100%. TEXT_SCALE multiplies every text element.

FIG_W_MM   <- 180                 # double-column width
TEXT_SCALE <- 1
PT <- function(x) x * TEXT_SCALE
SZ <- list(tick = 7, axis = 8, strip = 8, title = 8.5, legend = 7, annot = 6.5,
           venn_set = 7.5, venn_count = 6.5, tag = 12)
fig_w     <- FIG_W_MM / 25.4
GEOM_TEXT <- PT(SZ$annot) / .pt

theme_pub <- theme_classic(base_size = PT(SZ$axis)) +
  theme(plot.title    = element_text(face = "bold", hjust = 0.5, size = PT(SZ$title),
                                     margin = margin(b = 3)),
        plot.subtitle = element_text(hjust = 0.5, size = PT(SZ$annot), colour = "grey30"),
        axis.title    = element_text(size = PT(SZ$axis)),
        axis.text     = element_text(colour = "grey10", size = PT(SZ$tick)),
        strip.text    = element_text(face = "bold", size = PT(SZ$strip), margin = margin(2, 2, 2, 2)),
        strip.background = element_rect(fill = "grey95", colour = NA),
        plot.caption  = element_text(size = PT(SZ$annot), colour = "grey40", hjust = 0),
        legend.position = "bottom", legend.title = element_blank(),
        legend.text   = element_text(size = PT(SZ$legend)),
        legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(3, "pt"),
        plot.title.position = "panel",          # titles centred over the data
        plot.caption.position = "plot",
        plot.margin   = margin(4, 6, 4, 4))
title_theme <- theme(plot.title = element_text(face = "bold", hjust = 0.5, size = PT(SZ$title)))
tag_theme   <- theme(plot.tag = element_text(face = "bold", size = PT(SZ$tag)))

save_fig <- function(p, name, w, h) {
  ggsave(file.path(out_dir, paste0(name, ".pdf")), p, width = w, height = h,
         device = cairo_pdf, limitsize = FALSE)
  ggsave(file.path(out_dir, paste0(name, ".png")), p, width = w, height = h,
         dpi = 600, limitsize = FALSE, bg = "white")
}

# ---- 4. Read transcriptomics ----------------------------------------------

TX <- list()
for (ln in names(LINES)) {
  drs <- if (ln == "SVG-A") SVG_DRUGS else names(DRUGS)
  for (dg in drs) {
    f <- file.path(analysis_dir, tx_file(LINES[[ln]], dg))
    TX[[ln]][[dg]] <- read_csv(f, show_col_types = FALSE) %>%
      filter(!is.na(padj), !is.na(log2FoldChange)) %>%
      mutate(gene = na_if(trimws(gene_symbol), ""),
             sig  = padj < PADJ & abs(log2FoldChange) > LFC)
    message(sprintf("%-7s %-9s DETs %5d  DEGs %5d", ln, dg,
                    sum(TX[[ln]][[dg]]$sig),
                    n_distinct(na.omit(TX[[ln]][[dg]]$gene[TX[[ln]][[dg]]$sig]))))
  }
}

sig_genes <- function(ln, dg) unique(na.omit(TX[[ln]][[dg]]$gene[TX[[ln]][[dg]]$sig]))
sig_tx    <- function(ln, dg) unique(TX[[ln]][[dg]]$ENSEMBL_TRANSCRIPT[TX[[ln]][[dg]]$sig])
unique_sets <- function(sets) setNames(lapply(names(sets), function(x)
  setdiff(sets[[x]], unique(unlist(sets[names(sets) != x])))), names(sets))

# gene-level log2FC = most significant transcript per gene
# sig_only = TRUE: choose among significant transcripts (concordance, panel E)
gene_lfc <- function(ln, dg, sig_only = FALSE) {
  d <- TX[[ln]][[dg]]
  if (sig_only) d <- filter(d, sig)
  d %>% filter(!is.na(gene)) %>% arrange(padj, pvalue) %>%
    distinct(gene, .keep_all = TRUE) %>% select(gene, log2FoldChange, padj)
}


# ---- 5. Drug-unique transcriptomic features --------------------------------

UNIQ <- list(); counts <- list(); near_tab <- list()
for (ln in names(TX)) {
  drs <- names(TX[[ln]])
  G  <- setNames(lapply(drs, sig_genes, ln = ln), drs)
  Tr <- setNames(lapply(drs, sig_tx,    ln = ln), drs)
  uG <- unique_sets(G); uT <- unique_sets(Tr)
  UNIQ[[ln]] <- uG
  for (dg in drs) {
    near <- unique(unlist(lapply(setdiff(drs, dg), function(o)
      TX[[ln]][[o]] %>%
        filter(gene %in% uG[[dg]], pvalue < NEAR_P, abs(log2FoldChange) > NEAR_LFC) %>%
        pull(gene))))
    near_tab[[paste(ln, dg)]] <- tibble(line = ln, drug = dg, gene = uG[[dg]],
                                        near_threshold = uG[[dg]] %in% near)
    counts[[paste(ln, dg)]] <- tibble(line = ln, drug = dg,
                                      layer  = c("Transcripts", "Genes"),
                                      total  = c(length(Tr[[dg]]), length(G[[dg]])),
                                      unique = c(length(uT[[dg]]), length(uG[[dg]])))
  }
}
near_tab <- bind_rows(near_tab)


# ---- 6. Phosphoproteomics and metabolomics (KNS-42) ------------------------

# all regulated phosphopeptides
P_sets <- sapply(names(DRUGS), function(dg)
  unique(na.omit(trimws(as.character(
    read_excel(file.path(PHOS_DIR, phos_file(dg)))$gene)))), simplify = FALSE)
uP <- unique_sets(P_sets)
message(sprintf("Phosphoproteins: %d total, %d single-drug",
                length(unique(unlist(P_sets))), length(unlist(uP))))

MF <- read_excel(METAB_XLSX, sheet = METAB_SHEET)
METAB <- sapply(names(DRUGS), function(dg) {
  cc <- METAB_CONTRAST[[dg]]
  tibble(metabolite = trimws(as.character(MF$Name)),
         log2FC = -as.numeric(MF[[sprintf("Log2 Fold Change: (%s) / (%s)", cc[2], cc[1])]]),
         padj   =  as.numeric(MF[[sprintf("Adj. P-value: (%s) / (%s)",     cc[2], cc[1])]])) %>%
    filter(!is.na(metabolite), metabolite != "", !is.na(padj), !is.na(log2FC)) %>%
    mutate(sig = padj < PADJ & abs(log2FC) > LFC)
}, simplify = FALSE)
M_sets <- lapply(METAB, function(d) unique(d$metabolite[d$sig]))
uM <- unique_sets(M_sets)

counts <- bind_rows(
  bind_rows(counts),
  tibble(line = "KNS-42", drug = names(DRUGS), layer = "Phosphoproteins",
         total = lengths(P_sets), unique = lengths(uP)),
  tibble(line = "KNS-42", drug = names(DRUGS), layer = "Metabolites",
         total = lengths(M_sets), unique = lengths(uM))) %>%
  mutate(shared = total - unique)
print(counts, n = Inf)


# ---- 7. Test 1: cross-layer overlap of drug-unique features (KNS-42) -------

drs_K <- names(DRUGS)
uni_K <- Reduce(intersect, lapply(drs_K, function(d) unique(na.omit(TX[["KNS-42"]][[d]]$gene))))
uP_K  <- lapply(uP, intersect, uni_K)
genes_P <- unlist(uP_K[drs_K], use.names = FALSE)
lab_P   <- match(rep(drs_K, lengths(uP_K[drs_K])), drs_K)
in_uG   <- sapply(drs_K, function(d) genes_P %in% UNIQ[["KNS-42"]][[d]])  # gene x drug
idx     <- seq_along(genes_P)

obs_layer  <- sum(in_uG[cbind(idx, lab_P)])
null_layer <- replicate(N_PERM_LAYER, sum(in_uG[cbind(idx, sample(lab_P))]))
p_layer    <- mean(null_layer >= obs_layer)
layer_hits <- bind_rows(lapply(drs_K, function(d)
  tibble(drug = d, gene = intersect(UNIQ[["KNS-42"]][[d]], uP_K[[d]]))))
message(sprintf("Cross-layer: observed %d, expected %.2f, p = %.3f",
                obs_layer, mean(null_layer), p_layer))


# ---- 8. Test 2: cross-line reproducibility of drug-unique genes ------------

uni_both <- Reduce(intersect, unlist(lapply(TX, function(x)
  lapply(x, function(d) unique(na.omit(d$gene)))), recursive = FALSE))
line_res <- list(); line_null <- list()
for (dg in SVG_DRUGS) {
  k <- intersect(UNIQ[["KNS-42"]][[dg]], uni_both)
  s <- intersect(UNIQ[["SVG-A"]][[dg]],  uni_both)
  obs <- length(intersect(k, s))
  nul <- replicate(N_PERM_LINE, sum(sample(uni_both, length(k)) %in% s))
  line_null[[dg]] <- nul
  line_res[[dg]] <- tibble(
    drug = dg, KNS_unique = length(k), SVG_unique = length(s),
    observed = obs, expected = length(k) * length(s) / length(uni_both),
    p_perm = mean(nul >= obs),
    p_hyper = phyper(obs - 1, length(s), length(uni_both) - length(s),
                     length(k), lower.tail = FALSE),
    genes = paste(sort(intersect(k, s)), collapse = ";"))
}
line_res <- bind_rows(line_res)
print(select(line_res, -genes))


# ---- 9. Cross-line directional concordance (binomial) ----------------------

conc <- bind_rows(lapply(SVG_DRUGS, function(dg) {
  j <- inner_join(gene_lfc("KNS-42", dg, TRUE), gene_lfc("SVG-A", dg, TRUE),
                  by = "gene", suffix = c("_KNS", "_SVG"))
  n <- nrow(j); k <- sum(sign(j$log2FoldChange_KNS) == sign(j$log2FoldChange_SVG))
  tibble(drug = dg, shared_genes = n, concordant = k,
         pct_concordant = round(100 * k / n, 1),
         binom_p = if (n) binom.test(k, n, 0.5)$p.value else NA_real_)
}))
print(conc)


# ---- 10. Supplementary tables ----------------------------------------------

write_csv(counts,   file.path(out_dir, "S_drug_unique_counts.csv"))
write_csv(near_tab, file.path(out_dir, "S_drug_unique_genes_near_threshold.csv"))
write_csv(bind_rows(lapply(names(uP), function(d) tibble(drug = d, phosphoprotein = uP[[d]]))),
          file.path(out_dir, "S_drug_unique_phosphoproteins.csv"))
write_csv(bind_rows(lapply(names(DRUGS), function(d)
  METAB[[d]] %>% filter(sig, metabolite %in% uM[[d]]) %>% mutate(drug = d, .before = 1))),
  file.path(out_dir, "S_drug_unique_metabolites.csv"))
write_csv(tibble(test = "KNS-42 drug-unique DEGs x drug-unique phosphoproteins",
                 observed = obs_layer, expected = mean(null_layer),
                 p_perm = p_layer, n_perm = N_PERM_LAYER,
                 genes = paste(layer_hits$drug, layer_hits$gene, sep = ":", collapse = ";")),
          file.path(out_dir, "S_cross_layer_test.csv"))
write_csv(line_res, file.path(out_dir, "S_cross_line_test.csv"))
write_csv(conc,     file.path(out_dir, "S_crossline_concordance_binomial.csv"))


# ---- 11. Panel A: drug-unique vs shared counts -----------------------------

facet_lev <- c("KNS-42\nGenes", "KNS-42\nPhosphoproteins", "KNS-42\nMetabolites", "SVG-A\nGenes")
pd_A <- counts %>%
  filter(layer != "Transcripts") %>%
  mutate(facet = factor(paste(line, layer, sep = "\n"), levels = facet_lev),
         drug  = factor(drug, levels = names(DRUGS))) %>%
  pivot_longer(c(unique, shared), names_to = "class", values_to = "n") %>%
  mutate(class = factor(class, levels = c("shared", "unique"),
                        labels = c("Shared with another drug", "Drug-unique")))

p_A <- ggplot(pd_A, aes(drug, n, fill = drug, alpha = class)) +
  geom_col(colour = "grey20", linewidth = 0.2, width = 0.74) +
  facet_wrap(~ facet, nrow = 1, scales = "free") +
  scale_fill_manual(values = DRUG_COLS, guide = "none") +
  scale_alpha_manual(values = c(0.3, 1)) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.05))) +
  labs(title = "Drug-unique and shared features by cell line and omic layer",
       x = NULL, y = "Features (n)",
       caption = sprintf("Genes: %s. Phosphoproteins: regulated phosphopeptides (p < 0.05). Metabolites: %s. SVG-A morphine excluded.",
                         thr_lab, thr_lab)) +
  guides(alpha = guide_legend(override.aes = list(fill = "grey40"))) +
  theme_pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold", size = PT(SZ$annot), lineheight = 0.95),
        strip.clip = "off", panel.spacing.x = unit(10, "pt"))
save_fig(p_A, "Fig6A_unique_counts", fig_w, 2.1)


# ---- 12. Panel B: threshold dependence -------------------------------------

pd_B <- near_tab %>%
  count(line, drug, near_threshold) %>%
  group_by(line, drug) %>% mutate(pct = 100 * n / sum(n), tot = sum(n)) %>% ungroup() %>%
  mutate(drug = factor(drug, levels = names(DRUGS)),
         line = factor(line, levels = names(LINES)),
         class = factor(near_threshold, levels = c(TRUE, FALSE),
                        labels = c("Sub-threshold with another drug", "Robustly drug-unique")))

p_B <- ggplot(pd_B, aes(drug, pct, fill = drug, alpha = class)) +
  geom_col(colour = "grey20", linewidth = 0.2, width = 0.74) +
  geom_text(data = distinct(pd_B, line, drug, tot),
            aes(drug, 103, label = comma(tot)), inherit.aes = FALSE, size = GEOM_TEXT, vjust = 0) +
  facet_grid(~ line, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = DRUG_COLS, guide = "none") +
  scale_alpha_manual(values = c(0.3, 1)) +
  scale_y_continuous(limits = c(0, 116), breaks = seq(0, 100, 25),
                     expand = expansion(mult = c(0, 0))) +
  labs(title = "Threshold dependence of unique genes", x = NULL,
       y = "Genes (%)",
       caption = sprintf("Sub-threshold: p < %.2g and |log2FC| > %g for any other drug.\nNumbers = drug-unique genes.",
                         NEAR_P, NEAR_LFC)) +
  guides(alpha = guide_legend(nrow = 2, override.aes = list(fill = "grey40"))) +
  theme_pub + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_fig(p_B, "Fig6B_threshold_dependence", fig_w / 2, 2.4)


# ---- 13. Panel C: overlap of drug-unique features versus chance ------------

ci <- function(x) quantile(x, c(0.05, 0.95))     # one-sided tests
pd_C <- bind_rows(
  tibble(group = "Cross-layer", test = "Cross-layer: KNS-42\ngenes vs phosphoproteins",
         obs = obs_layer, exp = mean(null_layer),
         lo = ci(null_layer)[1], hi = ci(null_layer)[2], p = p_layer),
  bind_rows(lapply(SVG_DRUGS, function(dg) {
    x <- line_res[line_res$drug == dg, ]
    tibble(group = "Cross-line", test = paste0("Cross-line: ", dg),
           obs = x$observed, exp = x$expected,
           lo = ci(line_null[[dg]])[1], hi = ci(line_null[[dg]])[2], p = x$p_perm)
  }))) %>%
  mutate(label = sprintf("%s\n%d vs %.1f, p = %s", test, obs, exp, trimws(formatC(p, format = "g", digits = 2))),
         across(c(obs, lo, hi), ~ .x / exp, .names = "{.col}_r"),
         group = factor(group, levels = c("Cross-layer", "Cross-line")),
         label = factor(label, levels = rev(label)))
write_csv(pd_C %>% select(group, test, obs, exp, lo, hi, p), file.path(out_dir, "Fig6C_overlap_tests.csv"))

p_C <- ggplot(pd_C, aes(y = label)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey55", linewidth = 0.3) +
  geom_segment(aes(x = lo_r, xend = hi_r, yend = label), colour = COL_NS, linewidth = 2.2) +
  geom_point(aes(x = 1), shape = 21, fill = "white", colour = "grey20", size = 1.6, stroke = 0.4) +
  geom_point(aes(x = obs_r, colour = p < 0.05), size = 2) +
  scale_colour_manual(values = c(`TRUE` = COL_UP, `FALSE` = "grey20"), guide = "none") +
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  labs(title = "Drug-unique overlap versus chance",
       x = "Observed / expected overlap", y = NULL,
       caption = sprintf("Grey bar: 5th-95th percentile of %s / %s permutations; open circle: expected;\nfilled: observed (red, p < 0.05). Cross-line: KNS-42 vs SVG-A.",
                         comma(N_PERM_LAYER), comma(N_PERM_LINE))) +
  theme_pub +
  theme(axis.text.y = element_text(size = PT(SZ$annot), lineheight = 0.95),
        strip.text.y = element_blank(), strip.background = element_blank(),
        panel.spacing.y = unit(8, "pt"))
save_fig(p_C, "Fig6C_overlap_tests", fig_w / 2, 2.6)


# ---- 14. Supplementary: drug-unique genes across all drugs -----------------
# Genes run left to right, grouped by the drug they are unique to (equal
# width per group); rows are the drugs.

heat_data <- function(ln) {
  drs <- names(TX[[ln]])
  lf  <- bind_rows(lapply(drs, function(d) gene_lfc(ln, d) %>% mutate(drug = d)))
  bind_rows(lapply(drs, function(u) {
    ord <- lf %>% filter(drug == u, gene %in% UNIQ[[ln]][[u]]) %>%
      arrange(desc(log2FoldChange)) %>% pull(gene)
    expand_grid(unique_to = u, gene = ord, drug = drs) %>%
      left_join(lf, by = c("gene", "drug")) %>%
      mutate(rank = match(gene, ord),
             block = sprintf("%s (n = %s)", u, comma(length(ord))))
  })) %>%
    mutate(drug = factor(drug, levels = rev(names(DRUGS))),
           block = factor(block, levels = unique(block)),
           lfc = pmax(pmin(log2FoldChange, 5), -5))
}
heat_plot <- function(ln) {
  ggplot(heat_data(ln), aes(rank, drug, fill = lfc)) +
    geom_raster() +
    facet_grid(. ~ block, scales = "free_x") +
    scale_fill_gradient2(low = COL_DOWN, mid = "white", high = COL_UP, midpoint = 0,
                         limits = c(-5, 5), breaks = c(-5, 0, 5), na.value = COL_NS,
                         name = expression(log[2]~"fold-change")) +
    scale_x_continuous(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
    labs(title = ln, x = NULL, y = NULL) +
    theme_pub +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = PT(SZ$annot)),
          axis.line = element_blank(),
          strip.text = element_text(face = "bold", size = PT(SZ$annot)),
          panel.spacing.x = unit(3, "pt"),
          legend.title = element_text(size = PT(SZ$axis), vjust = 0.8),
          legend.key.height = unit(6, "pt"), legend.key.width = unit(40, "pt"))
}
make_heat <- function(caption) (heat_plot("KNS-42") / heat_plot("SVG-A")) +
  plot_layout(guides = "collect") +
  plot_annotation(title = "Response of drug-unique genes to every drug", caption = caption,
                  theme = title_theme + theme(plot.caption = element_text(size = PT(SZ$annot),
                                                                          colour = "grey40", hjust = 0))) &
  theme(legend.position = "bottom")
p_HEAT <- make_heat("Columns: drug-unique genes, ordered by log2FC for the drug they are unique to. Values capped at +/-5; grey = not quantified.")
if (SAVE_SUPPLEMENTARY) save_fig(p_HEAT, "FigS6A_unique_gene_heatmap", fig_w, 2.8)


# ---- 15. Supplementary: psilocin genes unique in both lines ----------------

psi_genes <- strsplit(line_res$genes[line_res$drug == "psilocin"], ";")[[1]]
pd_E <- inner_join(gene_lfc("KNS-42", "psilocin", TRUE),
                   gene_lfc("SVG-A",  "psilocin", TRUE),
                   by = "gene", suffix = c("_KNS", "_SVG")) %>%
  filter(gene %in% psi_genes) %>%
  mutate(dir = ifelse(sign(log2FoldChange_KNS) == sign(log2FoldChange_SVG),
                      "Concordant", "Discordant"))
n_dir <- table(pd_E$dir)
dir_lab <- setNames(sprintf("%s (n = %d)", names(n_dir), as.integer(n_dir)), names(n_dir))

p_E <- ggplot(pd_E, aes(log2FoldChange_KNS, log2FoldChange_SVG, colour = dir)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_point(size = 1.2, alpha = 0.85, stroke = 0) +
  geom_text_repel(aes(label = gene), size = GEOM_TEXT, fontface = "italic",
                  max.overlaps = Inf, seed = 1, show.legend = FALSE,
                  segment.size = 0.2, min.segment.length = 0.1,
                  box.padding = 0.12, point.padding = 0.05, force = 2) +
  scale_colour_manual(values = c(Concordant = COL_CONC, Discordant = COL_DISC),
                      labels = dir_lab) +
  guides(colour = guide_legend(override.aes = list(size = 2.5))) +
  labs(title = sprintf("Psilocin-unique genes in both lines (n = %d)", nrow(pd_E)),
       x = expression(log[2]~"fold-change"~"(KNS-42)"),
       y = expression(log[2]~"fold-change"~"(SVG-A)"),
       caption = sprintf("%s. Gene log2FC from the most significant transcript.", thr_lab)) +
  theme_pub
if (SAVE_SUPPLEMENTARY) save_fig(p_E, "FigS6B_psilocin_crossline", fig_w / 2, 3.0)


# ---- 16. Supplementary: psilocin-unique metabolites -----------------------

pd_F <- METAB$psilocin %>%
  filter(sig, metabolite %in% uM$psilocin) %>%
  group_by(metabolite) %>% slice_min(padj, n = 1, with_ties = FALSE) %>% ungroup() %>%
  mutate(label = ifelse(nchar(metabolite) > 35, paste0(substr(metabolite, 1, 33), "..."), metabolite),
         label = reorder(label, log2FC),
         cls = ifelse(log2FC > 0, "Up", "Down"))

p_F <- ggplot(pd_F, aes(log2FC, label, fill = cls)) +
  geom_col(colour = "grey20", linewidth = 0.2, width = 0.74) +
  geom_vline(xintercept = 0, colour = "grey20", linewidth = 0.3) +
  scale_fill_manual(values = c(Up = COL_UP, Down = COL_DOWN), breaks = c("Up", "Down")) +
  labs(title = "Psilocin-unique metabolites",
       x = expression(log[2]~"fold-change"), y = NULL, caption = thr_lab) +
  theme_pub + theme(axis.text.y = element_text(size = PT(SZ$annot)),
                    plot.margin = margin(4, 26, 4, 4))
if (SAVE_SUPPLEMENTARY) save_fig(p_F, "FigS6C_psilocin_metabolites", fig_w / 2, 3.0)


# ---- 17. Combined figure ---------------------------------------------------

# captions are dropped in the combined figure (they belong in the figure legend)
nocap <- function(p) p + labs(caption = NULL)
fig6 <- wrap_elements(full = nocap(p_A)) /
  (wrap_elements(full = nocap(p_B)) | wrap_elements(full = nocap(p_C))) +
  plot_layout(heights = c(2.0, 2.4)) +
  plot_annotation(tag_levels = list(c("A", "B", "C"))) &
  tag_theme
save_fig(fig6, "Fig6_combined", fig_w, 4.5)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
message("Figure 6 written to: ", out_dir)
