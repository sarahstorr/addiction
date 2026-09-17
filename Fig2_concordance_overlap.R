

# ---- 1. Settings ----------------------------------------------------------

analysis_dir <- "C:/Users/mrzsjs1/Desktop/Addiction/analysis"
out_dir      <- file.path(analysis_dir, "figures", "Fig2")

LINES <- c("KNS-42" = "KNS", "SVG-A" = "SVG")       # cell line -> filename prefix
DRUGS <- c(ethanol = "ethanol", THC = "thc", cocaine = "cocaine",
           morphine = "morphine", psilocin = "psilocin")

EXCLUDE <- list(c("SVG-A", "morphine"))

## KNS-42 ethanol and psilocin use the matched-control reruns (*_CORRECTED files)
USE_CORRECTED <- TRUE

# column aliases (matched case-insensitively)
ID_ALIASES   <- c("ensembl_transcript_ID_version","ENSEMBL_TRANSCRIPT","ensembl_transcript",
                  "transcript_id","target_id","...1","X")
LFC_ALIASES  <- c("log2FoldChange","log2FC","log2fc","logFC","log2_fold_change")
PADJ_ALIASES <- c("padj","padj_BH","adj.P.Val","FDR","p_adj","qvalue")
GENE_ALIASES <- c("gene_symbol","external_gene_name","gene_name","symbol")

# significance rule - SINGLE SOURCE OF TRUTH
PADJ <- 0.05; LFC <- 1
STRIP_VERSION <- FALSE

# palette (unchanged from the original scripts)
COL_CONC   <- "#1D9E75"; COL_DISC <- "#E24B4A"
VENN_LOW   <- "#F7FBFF"; VENN_HIGH <- "#7F77DD"
SHOW_CONC_LEGEND <- TRUE          # concordant/discordant key under panels A and B



# ---- 2. Packages ----------------------------------------------------------

pkgs <- c("readr","dplyr","tidyr","ggplot2","scales","patchwork","ggVennDiagram")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing)
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
  library(scales); library(patchwork); library(ggVennDiagram)
})
# guard against masking by Bioconductor packages still attached
select <- dplyr::select; filter <- dplyr::filter
rename <- dplyr::rename; mutate <- dplyr::mutate

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
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

row_title <- function(txt) plot_annotation(title = txt, theme = title_theme)

# ---- 4. Read transcriptomic results --------------------------------------

pick <- function(aliases, nm, label) {
  hit <- nm[tolower(nm) %in% tolower(aliases)]
  if (!length(hit)) stop("[", label, "] no column among: ", paste(aliases, collapse = ", "),
                         "\n  file has: ", paste(nm, collapse = ", "))
  hit[1]
}
find_file <- function(line_key, drug_key) {
  f <- list.files(analysis_dir, pattern = "\\.csv$", full.names = TRUE)
  f <- f[!grepl("noGeneSymbol", basename(f), ignore.case = TRUE)]
  hit <- f[grepl(paste0("(^|[^A-Za-z])", line_key), basename(f), ignore.case = TRUE) &
             grepl(drug_key, basename(f), ignore.case = TRUE)]
  if (!length(hit)) return(NA_character_)
  corr <- grepl("CORRECTED", basename(hit), ignore.case = TRUE)
  ## prefer the matched-control rerun where one exists, else the standard file
  if (USE_CORRECTED && any(corr)) hit[corr][1] else hit[!corr][1]
}
read_full <- function(line, drug) {
  fp <- find_file(LINES[[line]], DRUGS[[drug]])
  if (is.na(fp)) { warning("missing file: ", line, " / ", drug); return(NULL) }
  df <- readr::read_csv(fp, show_col_types = FALSE); nm <- names(df); lab <- paste(line, drug)
  cid <- pick(ID_ALIASES, nm, lab); clf <- pick(LFC_ALIASES, nm, lab); cpj <- pick(PADJ_ALIASES, nm, lab)
  cgene <- GENE_ALIASES[tolower(GENE_ALIASES) %in% tolower(nm)]
  cgene <- if (length(cgene)) cgene[1] else NA_character_
  d <- df %>% transmute(id = as.character(.data[[cid]]),
                        lfc = as.numeric(.data[[clf]]),
                        padj = as.numeric(.data[[cpj]]),
                        gene = if (!is.na(cgene)) toupper(trimws(as.character(.data[[cgene]]))) else NA_character_) %>%
    filter(!is.na(id), !is.na(lfc), !is.na(padj))
  if (STRIP_VERSION) d$id <- sub("\\..*$", "", d$id)
  d$gene <- dplyr::na_if(d$gene, "")
  d$drug <- drug; d$line <- line
  d$sig <- d$padj < PADJ & abs(d$lfc) > LFC
  message(sprintf("  %-14s %5d DETs, %5d DEGs  <- %s", lab, sum(d$sig),
                  dplyr::n_distinct(d$gene[d$sig & !is.na(d$gene)]), basename(fp)))
  d
}

## only the significant rows are needed for Figure 2
SIG <- list()
for (ln in names(LINES)) for (dg in names(DRUGS)) {
  if (any(vapply(EXCLUDE, function(e) e[1] == ln && e[2] == dg, logical(1)))) {
    message(sprintf("  %-14s excluded", paste(ln, dg))); next
  }
  d <- read_full(ln, dg)
  if (!is.null(d)) SIG[[paste(ln, dg)]] <- dplyr::filter(d, sig)
}
SIG <- bind_rows(SIG)
if (!nrow(SIG)) stop("No CSVs read - check analysis_dir and filename keywords.")

## one log2FC per DET (transcript) or per DEG (most significant DET of that gene)
level_values <- function(ln, dg, level) {
  d <- SIG %>% filter(line == ln, drug == dg)
  if (level == "transcript") return(d %>% transmute(key = id, lfc))
  d %>% filter(!is.na(gene)) %>% group_by(gene) %>%
    slice_min(padj, n = 1, with_ties = FALSE) %>% ungroup() %>%
    transmute(key = gene, lfc)
}


# ---- 5. Panels A and B: cross-line concordance ---------------------------

## a single key row, drawn directly. Legends collected by patchwork can appear
## more than once when the panels differ, so the key is built here instead.
legend_row <- function(labels, colours) {
  n <- length(labels)
  d <- data.frame(x = (seq_len(n) - 1) * 1.6, lab = labels, col = colours)
  ggplot(d, aes(x, 1)) +
    geom_point(aes(colour = lab), size = 2.5, show.legend = FALSE) +
    geom_text(aes(label = lab), hjust = 0, nudge_x = 0.14,
              size = PT(SZ$legend) / .pt, colour = "grey10") +
    scale_colour_manual(values = setNames(colours, labels)) +
    ## symmetric limits keep the keys clustered at the centre of the row
    scale_x_continuous(limits = c(-1.5, (n - 1) * 1.6 + 1.5)) +
    coord_cartesian(clip = "off") +
    theme_void() + theme(plot.margin = margin(0, 0, 0, 0))
}

make_conc <- function(mg, ttl) {
  axm <- max(abs(c(mg$lfc_K, mg$lfc_S)))
  ggplot(mg, aes(lfc_K, lfc_S, colour = dir)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
    geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.3) +
    geom_point(size = 1.2, alpha = 0.85, stroke = 0) +
    scale_colour_manual(values = c(Concordant = COL_CONC, Discordant = COL_DISC),
                        drop = FALSE, guide = "none") +
    coord_equal(xlim = c(-axm, axm), ylim = c(-axm, axm)) +
    labs(title = ttl,
         x = expression(log[2]*"FC, KNS-42"), y = expression(log[2]*"FC, SVG-A")) +
    theme_pub
}

conc_row <- function(level = c("transcript","gene")) {
  level <- match.arg(level)
  ps <- list(); stats <- list()
  for (dg in names(DRUGS)) {
    k <- level_values("KNS-42", dg, level); s <- level_values("SVG-A", dg, level)
    if (!nrow(k) || !nrow(s)) next                       # e.g. SVG-A morphine
    mg <- inner_join(k %>% rename(lfc_K = lfc), s %>% rename(lfc_S = lfc), by = "key") %>%
      mutate(dir = factor(ifelse(sign(lfc_K) == sign(lfc_S), "Concordant", "Discordant"),
                          levels = c("Concordant","Discordant")))
    if (!nrow(mg)) next
    n_t <- nrow(mg); n_c <- sum(mg$dir == "Concordant")
    pr  <- if (n_t > 2) cor(mg$lfc_K, mg$lfc_S) else NA_real_
    ps[[dg]] <- make_conc(mg, dg)
    stats[[dg]] <- data.frame(drug = dg, level = level, shared = n_t, concordant = n_c,
                              pct_concordant = round(100 * n_c / n_t, 1),
                              pearson_r = round(pr, 3))
    write_csv(mg %>% rename(!!level := key),
              file.path(out_dir, sprintf("Fig2_%s_%s_shared.csv", level, dg)))
  }
  p <- wrap_plots(ps, nrow = 1)
  if (SHOW_CONC_LEGEND)
    p <- p / legend_row(c("Concordant", "Discordant"), c(COL_CONC, COL_DISC)) +
      plot_layout(heights = c(1, 0.1))
  list(plot = p, stats = bind_rows(stats))
}

c_det <- conc_row("transcript"); c_deg <- conc_row("gene")
conc_stats <- bind_rows(c_det$stats, c_deg$stats)
print(conc_stats, row.names = FALSE)
write_csv(conc_stats, file.path(out_dir, "Fig2AB_concordance_stats.csv"))

p_A <- c_det$plot + row_title("Shared-transcript overlap & directional concordance")
p_B <- c_deg$plot + row_title("Shared-gene overlap & directional concordance")
save_fig(p_A, "Fig2A_transcript_concordance", fig_w, 2.1)
save_fig(p_B, "Fig2B_gene_concordance",       fig_w, 2.1)


# ---- 6. Panels C and D: within-line overlap across drugs -----------------

set_list <- function(ln, level) {
  out <- list()
  for (dg in names(DRUGS)) {
    v <- level_values(ln, dg, level)$key
    if (length(v)) out[[dg]] <- unique(v)
  }
  out
}

venn_pair <- function(level = c("transcript","gene")) {
  level <- match.arg(level)
  unit_lab <- if (level == "transcript") "DETs" else "DEGs"
  ps <- list(); counts <- list()
  for (ln in names(LINES)) {
    sets <- set_list(ln, level)
    if (length(sets) < 2) next
    ps[[ln]] <- make_venn(sets, sprintf("%s: shared %s", ln, unit_lab),
                          unit_lab)
    ## private (drug-only) counts, used to check against the manuscript figure
    counts[[ln]] <- data.frame(line = ln, level = level, drug = names(sets),
                               total = lengths(sets),
                               private = vapply(names(sets), function(d)
                                 length(setdiff(sets[[d]], unlist(sets[names(sets) != d]))),
                                 integer(1)))
  }
  list(plot = wrap_plots(ps, nrow = 1), counts = bind_rows(counts))
}

v_det <- venn_pair("transcript"); v_deg <- venn_pair("gene")
venn_counts <- bind_rows(v_det$counts, v_deg$counts)
print(venn_counts, row.names = FALSE)
write_csv(venn_counts, file.path(out_dir, "Fig2CD_venn_counts.csv"))

save_fig(v_det$plot, "Fig2C_DET_venn", fig_w, 2.5)
save_fig(v_deg$plot, "Fig2D_DEG_venn", fig_w, 2.5)


# ---- 7. Assemble Figure 2 ------------------------------------------------

fig2 <- wrap_elements(full = p_A) /
  wrap_elements(full = p_B) /
  wrap_elements(full = v_det$plot) /
  wrap_elements(full = v_deg$plot) +
  plot_layout(heights = c(2.1, 2.1, 2.5, 2.5)) +
  plot_annotation(tag_levels = list(c("A", "B", "C", "D"))) &
  tag_theme
save_fig(fig2, "Fig2_combined", fig_w, 9.2)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
message("Figure 2 written to: ", out_dir)
