
# ---- 1. Settings ----------------------------------------------------------

analysis_dir <- "C:/Users/mrzsjs1/Desktop/Addiction/analysis"
out_dir      <- file.path(analysis_dir, "figures", "Fig4")

METAB_XLSX  <- "C:/Users/mrzsjs1/Desktop/Addiction/analysis/metabolomics final.xlsx"
METAB_SHEET <- "Filtered"
DET_DIR     <- analysis_dir                         # annotated KNS DESeq2 CSVs
PHOS_DIR    <- "C:/Users/mrzsjs1/Downloads"          # *_phosphoproteomics_significant.xlsx

DRUGS <- c(Ethanol = "ethanol", THC = "thc", Cocaine = "cocaine",
           Morphine = "morphine", Psilocin = "psilocin")      # display order
## metabolomics contrasts: drug group vs its control (THC vs ethanol vehicle)
CONTRAST <- list(Ethanol = c("E","UC"), THC = c("T","EC"), Cocaine = c("C","UC"),
                 Morphine = c("M","UC"), Psilocin = c("P","UC"))

PADJ <- 0.05; LFC <- 1

## panel C
LINE_KEY       <- "KNS"
USE_CORRECTED  <- TRUE       # KNS ethanol/psilocin matched-control reruns
PHOS_USE_RAW_P <- TRUE       # phosphopeptides: raw p (few survive correction)
COLLAPSE <- "most_sig"       # one value per gene: "most_sig" (smallest p) or "mean"
GENESET  <- "either"         # genes significant in "either" layer, "both", or "all"
SHOW_YX  <- TRUE             # dashed y = x line

# palette
COL_UP <- "#E24B4A"; COL_DOWN <- "#378ADD"; COL_NS <- "#B4B2A9"
VENN_LOW <- "#F7FBFF"; VENN_HIGH <- "#7F77DD"
COL_CONC_UP <- "#1D9E75"; COL_CONC_DOWN <- "#378ADD"; COL_DISC <- "#E69F00"


# ---- 2. Packages ----------------------------------------------------------

pkgs <- c("readr","readxl","dplyr","tidyr","ggplot2","scales","patchwork","ggVennDiagram")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing)
suppressPackageStartupMessages({
  library(readr); library(readxl); library(dplyr); library(tidyr); library(ggplot2)
  library(scales); library(patchwork); library(ggVennDiagram)
})
select <- dplyr::select; filter <- dplyr::filter; rename <- dplyr::rename
mutate <- dplyr::mutate; slice <- dplyr::slice; arrange <- dplyr::arrange
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
if (!file.exists(METAB_XLSX)) stop("Metabolomics file not found: ", METAB_XLSX)

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
# ---- Venn helper ---------------------------------------------------------------
# Set names are placed just outside the outlines, pointing away from the centre,
# so they never overlap the shapes. Falls back to ggVennDiagram's own placement
# if the installed version lacks process_data()/venn_setlabel().

venn_set_labels <- function(sets) {
  tryCatch({
    vd <- ggVennDiagram::process_data(ggVennDiagram::Venn(sets))
    ed <- ggVennDiagram::venn_setedge(vd)
    lb <- ggVennDiagram::venn_setlabel(vd)
    cx <- mean(range(ed$X)); cy <- mean(range(ed$Y))
    span <- max(diff(range(ed$X)), diff(range(ed$Y)))
    ae <- atan2(ed$Y - cy, ed$X - cx)
    re <- sqrt((ed$X - cx)^2 + (ed$Y - cy)^2)
    do.call(rbind, lapply(seq_len(nrow(lb)), function(i) {
      a    <- atan2(lb$Y[i] - cy, lb$X[i] - cx)
      near <- abs(atan2(sin(ae - a), cos(ae - a))) < 0.2
      R    <- max(re[near]) + 0.04 * span
      data.frame(name = as.character(lb$name[i]),
                 x = cx + R * cos(a), y = cy + R * sin(a),
                 hjust = if (cos(a) > 0.3) 0 else if (cos(a) < -0.3) 1 else 0.5,
                 vjust = if (sin(a) > 0.3) 0 else if (sin(a) < -0.3) 1 else 0.5)
    }))
  }, error = function(e) NULL)
}

make_venn <- function(sets, ttl, fill_name) {
  lab <- venn_set_labels(sets)
  p <- ggVennDiagram(sets, label = "count", label_alpha = 0, edge_size = 0.35,
                     set_size   = if (is.null(lab)) PT(SZ$venn_set) / .pt else 0,
                     label_size = PT(SZ$venn_count) / .pt) +
    scale_fill_gradient(low = VENN_LOW, high = VENN_HIGH, name = fill_name, labels = comma)
  if (!is.null(lab))
    p <- p + geom_text(data = lab, aes(x = x, y = y, label = name, hjust = hjust, vjust = vjust),
                       inherit.aes = FALSE, size = PT(SZ$venn_set) / .pt)
  p +
    scale_x_continuous(expand = expansion(mult = 0.22)) +
    scale_y_continuous(expand = expansion(mult = 0.10)) +
    coord_equal(clip = "off") +
    labs(title = ttl) +
    theme(plot.title   = element_text(face = "bold", hjust = 0.5, size = PT(SZ$title),
                                      margin = margin(b = 3)),
          plot.title.position = "panel",
          legend.position = "right",
          legend.title = element_text(size = PT(SZ$legend)),
          legend.text  = element_text(size = PT(SZ$legend)),
          legend.key.height = unit(13, "pt"), legend.key.width = unit(6, "pt"),
          legend.box.spacing = unit(12, "pt"),
          plot.margin = margin(4, 4, 4, 4))
}

# ---- 4. Read metabolomics --------------------------------------------------

## Compound Discoverer reports ratios as control / drug; log2FC is negated to
## give drug vs control (e.g. benzoylecgonine is strongly positive in cocaine)
MF <- read_excel(METAB_XLSX, sheet = METAB_SHEET)
METAB <- bind_rows(lapply(names(CONTRAST), function(d) {
  cc <- CONTRAST[[d]]
  tibble(Name   = as.character(MF$Name),
         log2FC = -as.numeric(MF[[sprintf("Log2 Fold Change: (%s) / (%s)", cc[2], cc[1])]]),
         adj_p  =  as.numeric(MF[[sprintf("Adj. P-value: (%s) / (%s)",     cc[2], cc[1])]]),
         drug   = d)
})) %>% filter(!is.na(Name), !is.na(log2FC), !is.na(adj_p)) %>%
  mutate(drug = factor(drug, levels = names(DRUGS)),
         sig  = adj_p < PADJ & abs(log2FC) > LFC,
         cls  = factor(ifelse(!sig, "ns", ifelse(log2FC > 0, "up", "down")),
                       levels = c("ns","down","up")))

mcounts <- METAB %>% filter(sig) %>% group_by(drug) %>%
  summarise(n_up = sum(log2FC > 0), n_dn = sum(log2FC < 0), .groups = "drop")
print(as.data.frame(mcounts), row.names = FALSE)
write_csv(mcounts, file.path(out_dir, "Fig4A_metabolite_counts.csv"))


# ---- 5. Panel A: volcano plots --------------------------------------------

mx   <- c(floor(min(METAB$log2FC)), ceiling(max(METAB$log2FC)))
ydat <- max(-log10(METAB$adj_p))
my   <- ceiling(ydat + 2)

p_A <- ggplot(METAB, aes(log2FC, -log10(adj_p))) +
  geom_hline(yintercept = -log10(PADJ), linetype = "dashed", colour = "grey55", linewidth = 0.3) +
  annotate("segment", x = c(-LFC, LFC), xend = c(-LFC, LFC), y = 0, yend = ydat + 0.5,
           linetype = "dashed", colour = "grey55", linewidth = 0.3) +
  geom_point(aes(colour = cls), size = 0.7, alpha = 0.7, stroke = 0) +
  geom_text(data = mcounts, aes(x = mx[1], y = my, label = paste(n_dn, "down")),
            inherit.aes = FALSE, hjust = 0, vjust = 1, size = GEOM_TEXT,
            colour = COL_DOWN, fontface = "bold") +
  geom_text(data = mcounts, aes(x = mx[2], y = my, label = paste(n_up, "up")),
            inherit.aes = FALSE, hjust = 1, vjust = 1, size = GEOM_TEXT,
            colour = COL_UP, fontface = "bold") +
  scale_colour_manual(values = c(ns = COL_NS, down = COL_DOWN, up = COL_UP), guide = "none") +
  scale_x_continuous(breaks = seq(-10, 20, 5)) +
  scale_y_continuous(breaks = seq(0, 30, 5), expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(xlim = mx, ylim = c(0, my)) +
  facet_wrap(~ drug, nrow = 1) +
  labs(title = "Metabolomic response by drug (KNS-42)",
       x = expression(log[2]~"fold-change"), y = expression(-log[10]~adjusted~italic(p)),
       caption = sprintf("Adj. p < %g and |log2FC| > %g; grey, not significant. THC is compared with the ethanol vehicle.",
                         PADJ, LFC)) +
  theme_pub
save_fig(p_A, "Fig4A_metabolomic_volcanoes", fig_w, 1.9)


# ---- 6. Panel B: metabolite features shared across drugs -------------------

metab_sets <- lapply(split(METAB$Name[METAB$sig], METAB$drug[METAB$sig]), unique)
metab_sets <- metab_sets[lengths(metab_sets) > 0]
p_B <- make_venn(metab_sets, "Shared metabolite features", "features")
venn_counts <- data.frame(drug = names(metab_sets), total = lengths(metab_sets),
  private = vapply(names(metab_sets), function(d)
    length(setdiff(metab_sets[[d]], unlist(metab_sets[names(metab_sets) != d]))), integer(1)),
  in_all = length(Reduce(intersect, metab_sets)))
print(venn_counts, row.names = FALSE)
write_csv(venn_counts, file.path(out_dir, "Fig4B_venn_counts.csv"))
save_fig(p_B, "Fig4B_metabolite_venn", 3.6, 2.9)


# ---- 7. Panel C: transcriptome vs phosphoproteome concordance -------------

pick <- function(aliases, nm, field, required = TRUE) {
  hit <- nm[tolower(nm) %in% tolower(aliases)]
  if (!length(hit)) {
    if (required) stop("No ", field, " column. Looked for: ", paste(aliases, collapse = ", "))
    return(NA_character_)
  }
  hit[1]
}
find_det <- function(drug_key) {
  f <- list.files(DET_DIR, pattern = "\\.csv$", full.names = TRUE)
  f <- f[!grepl("noGeneSymbol", basename(f), ignore.case = TRUE)]
  hit <- f[grepl(paste0("(^|[^A-Za-z])", LINE_KEY), basename(f), ignore.case = TRUE) &
             grepl(drug_key, basename(f), ignore.case = TRUE)]
  if (!length(hit)) return(NA_character_)
  corr <- grepl("CORRECTED", basename(hit), ignore.case = TRUE)
  if (USE_CORRECTED && any(corr)) hit[corr][1] else hit[!corr][1]
}
find_phos <- function(drug_key) {
  f <- list.files(PHOS_DIR, pattern = "\\.xlsx?$", full.names = TRUE)
  hit <- f[grepl(drug_key, basename(f), ignore.case = TRUE) & grepl("phospho", basename(f), ignore.case = TRUE)]
  if (length(hit)) hit[1] else NA_character_
}
collapse_layer <- function(df, val, pcol) {
  if (COLLAPSE == "mean") {
    df %>% group_by(gene) %>%
      summarise(v = mean(.data[[val]], na.rm = TRUE),
                p = suppressWarnings(min(.data[[pcol]], na.rm = TRUE)), .groups = "drop")
  } else {
    df %>% arrange(.data[[pcol]]) %>% group_by(gene) %>% slice(1) %>% ungroup() %>%
      transmute(gene, v = .data[[val]], p = .data[[pcol]])
  }
}
build_drug <- function(drug_label, drug_key) {
  det_fp <- find_det(drug_key); phos_fp <- find_phos(drug_key)
  if (is.na(det_fp) || is.na(phos_fp)) { message("  ", drug_label, ": missing file - skipped"); return(NULL) }
  dr <- read_csv(det_fp, show_col_types = FALSE); nm <- names(dr)
  dg <- pick(c("gene_symbol","external_gene_name","gene_name","symbol"), nm, "DET gene")
  dl <- pick(c("log2FoldChange","log2FC"), nm, "DET log2FC")
  dp <- pick(c("padj","FDR","adj.P.Val"), nm, "DET padj")
  det <- dr %>% transmute(gene = toupper(trimws(as.character(.data[[dg]]))),
                          t_lfc = as.numeric(.data[[dl]]), t_p = as.numeric(.data[[dp]])) %>%
    filter(!is.na(gene), gene != "", !is.na(t_lfc), !is.na(t_p))
  det_c <- collapse_layer(det, "t_lfc", "t_p") %>% rename(t_lfc = v, t_p = p)

  pr <- read_excel(phos_fp); pm <- names(pr)
  pg   <- pick(c("gene","gene_name","external_gene_name","symbol"), pm, "phospho gene")
  pl   <- pick(c("log2_fold_change","log2FC","log2FoldChange","logFC"), pm, "phospho log2FC")
  praw <- pick(c("p_value","pvalue","P.Value"), pm, "phospho raw p", required = FALSE)
  padj <- pick(c("adj_p_value","padj","adj.P.Val","FDR"), pm, "phospho adj p", required = FALSE)
  psig <- if (PHOS_USE_RAW_P && !is.na(praw)) praw else padj
  phos <- pr %>% transmute(gene = toupper(trimws(sub("[;,].*$", "", as.character(.data[[pg]])))),
                           p_lfc = as.numeric(.data[[pl]]),
                           p_p = if (!is.na(psig)) as.numeric(.data[[psig]]) else NA_real_) %>%
    filter(!is.na(gene), gene != "", !is.na(p_lfc))
  phos_c <- collapse_layer(phos, "p_lfc", "p_p") %>% rename(p_lfc = v, p_p = p)

  mg <- inner_join(det_c, phos_c, by = "gene")
  det_sig  <- !is.na(mg$t_p) & mg$t_p < PADJ & abs(mg$t_lfc) > LFC
  phos_sig <- if (all(is.na(mg$p_p))) abs(mg$p_lfc) > LFC else (mg$p_p < PADJ & abs(mg$p_lfc) > LFC)
  keep <- switch(GENESET, both = det_sig & phos_sig, either = det_sig | phos_sig, all = TRUE)
  mg <- mg[keep, , drop = FALSE]
  if (!nrow(mg)) return(NULL)
  mg$class <- with(mg, ifelse(sign(t_lfc) == sign(p_lfc),
                              ifelse(t_lfc > 0, "Concordant up", "Concordant down"), "Discordant"))
  mg$drug <- drug_label
  message(sprintf("  %-9s %3d genes <- %s", drug_label, nrow(mg), basename(det_fp)))
  mg
}

dat <- lapply(names(DRUGS), function(d) build_drug(d, DRUGS[[d]])); names(dat) <- names(DRUGS)
dat <- dat[!vapply(dat, is.null, logical(1))]
if (!length(dat)) stop("No drugs produced data for panel C - check DET_DIR / PHOS_DIR.")

col_map <- c("Concordant up" = COL_CONC_UP, "Concordant down" = COL_CONC_DOWN, "Discordant" = COL_DISC)
ax_all  <- max(vapply(dat, function(m) max(abs(c(m$t_lfc, m$p_lfc))), numeric(1)))

one_panel <- function(mg, first = TRUE) {
  rr <- if (nrow(mg) > 2) cor(mg$t_lfc, mg$p_lfc) else NA_real_
  ggplot(mg, aes(t_lfc, p_lfc, colour = factor(class, levels = names(col_map)))) +
    geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.4) +
    geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.4) +
    { if (SHOW_YX) geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                               colour = "grey55", linewidth = 0.4) } +
    geom_point(size = 1.1, alpha = 0.8, stroke = 0) +
    scale_colour_manual(values = col_map, drop = FALSE) +
    guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
    scale_x_continuous(breaks = c(-10, 0, 10)) +
    scale_y_continuous(breaks = c(-10, -5, 0, 5, 10)) +
    coord_equal(xlim = c(-ax_all, ax_all), ylim = c(-ax_all, ax_all)) +
    labs(title = mg$drug[1],
         subtitle = sprintf("%d genes\nr = %s", nrow(mg),
                            ifelse(is.na(rr), "NA", formatC(rr, format = "f", digits = 2))),
         x = expression(RNA~log[2]*FC),
         y = if (first) expression(Phospho~log[2]*FC) else NULL) +
    theme_pub +
    theme(plot.subtitle = element_text(lineheight = 0.95, hjust = 0.5, size = PT(SZ$annot), colour = "grey30"))
}
p_C <- wrap_plots(lapply(seq_along(dat), function(i) one_panel(dat[[i]], first = i == 1)), nrow = 1) +
  plot_layout(guides = "collect") +
  plot_annotation(title = "Transcriptomic vs phosphoproteomic concordance in KNS-42 by drug",
                  theme = title_theme) &
  theme(legend.position = "bottom")
save_fig(p_C, "Fig4C_transcript_phospho_concordance", fig_w, 2.3)

conc <- bind_rows(lapply(dat, function(m) data.frame(
  drug = m$drug[1], genes = nrow(m),
  conc_up = sum(m$class == "Concordant up"), conc_down = sum(m$class == "Concordant down"),
  discordant = sum(m$class == "Discordant"), pearson_r = round(cor(m$t_lfc, m$p_lfc), 3))))
print(conc, row.names = FALSE)
write_csv(conc, file.path(out_dir, "Fig4C_concordance_summary.csv"))


# ---- 8. Assemble Figure 4 --------------------------------------------------

row_B <- plot_spacer() + wrap_elements(full = p_B) + plot_spacer() +
  plot_layout(widths = c(1, 1.6, 1))                  # Venn centred on the page
fig4 <- wrap_elements(full = p_A) /
  row_B /
  wrap_elements(full = p_C) +
  plot_layout(heights = c(1.9, 2.9, 2.3)) +
  plot_annotation(tag_levels = list(c("A", "B", "C"))) &
  tag_theme
save_fig(fig4, "Fig4_combined", fig_w, 7.2)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
message("Figure 4 written to: ", out_dir)
