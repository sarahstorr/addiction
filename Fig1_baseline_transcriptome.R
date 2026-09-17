
# ---- 1. Settings ----------------------------------------------------------

analysis_dir <- "C:/Users/mrzsjs1/Desktop/Addiction/analysis"
out_dir      <- file.path(analysis_dir, "figures", "Fig1")

## Panel B: the GFAP immunofluorescence image.

IF_IMAGE <- file.path(analysis_dir, "GFAP.png")

LINES <- c("KNS-42" = "KNS", "SVG-A" = "SVG")       # cell line -> filename prefix
DRUGS <- c(ethanol = "ethanol", THC = "thc", cocaine = "cocaine",
           morphine = "morphine", psilocin = "psilocin")
VENN_DRUGS <- c("ethanol", "THC", "cocaine", "psilocin")   # SVG-A morphine excluded


EXCLUDE <- list(c("SVG-A", "morphine"))

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
Y_CAP <- 60                       # volcano -log10(padj) display cap (NA = none)

# palette (unchanged from the original scripts)
COL_UP <- "#E24B4A"; COL_DOWN <- "#378ADD"; COL_NS <- "#B4B2A9"
VENN2_COLS <- c("#1D9E75", "#7F77DD")               # KNS-42 teal / SVG-A purple
PROT_FILL  <- "#1D9E75"
CELL_LINE  <- "KNS-42"
SHOW_TIERS <- FALSE               # expression tier brackets above panel A

# font: original figures used base 15 / title 17 / axis text 13


# ---- 2. Packages ----------------------------------------------------------

pkgs <- c("readr","dplyr","tidyr","ggplot2","scales","patchwork","png","jpeg","grid")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing)
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
  library(scales); library(patchwork); library(grid)
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

# ---- 4. Panel A: baseline astrocyte marker proteome ----------------------

prot <- data.frame(
  gene_name = c("VIM","LGALS3BP","GFAP","PEA15","GAP43","SLC1A3",
                "AQP4","GJA1","GSS","SLC1A2","NOTCH1"),
  log10_abundance = c(8.641224791,6.244483122,6.932189985,6.218629496,6.100859465,
                      5.674796821,5.607878246,5.550693136,5.427245577,5.261796419,4.681234142),
  expression_level = c("Very High","Very High","Very High","Very High","Very High",
                       "High","High","High","High","High","Medium"),
  stringsAsFactors = FALSE) %>%
  arrange(desc(log10_abundance)) %>%
  mutate(gene_name = factor(gene_name, levels = gene_name),
         expression_level = factor(expression_level, levels = c("Very High","High","Medium")))

y_bot <- 4; y_top <- ceiling(max(prot$log10_abundance))
tier_layer <- list()
if (SHOW_TIERS) {
  tiers <- prot %>% mutate(idx = as.integer(gene_name)) %>%
    group_by(expression_level) %>%
    summarise(x1 = min(idx) - 0.4, x2 = max(idx) + 0.4, .groups = "drop") %>%
    mutate(xmid = (x1 + x2) / 2)
  ybr <- y_top + 0.22
  tier_layer <- list(
    geom_segment(data = tiers, aes(x = x1, xend = x2, y = ybr, yend = ybr),
                 inherit.aes = FALSE, colour = "grey60", linewidth = 0.35),
    geom_text(data = tiers, aes(x = xmid, y = ybr + 0.12, label = expression_level),
              inherit.aes = FALSE, size = GEOM_TEXT, colour = "grey45", fontface = "italic"))
  y_top <- y_top + 0.55
}

p_prot <- ggplot(prot, aes(gene_name, log10_abundance)) +
  geom_col(fill = PROT_FILL, colour = "grey20", linewidth = 0.25, width = 0.74) +
  tier_layer +
  scale_y_continuous(breaks = seq(y_bot, ceiling(max(prot$log10_abundance)), 1),
                     labels = function(b) parse(text = paste0("10^", b)),
                     expand = expansion(mult = c(0, 0.02)),
                     limits = c(y_bot, y_top), oob = scales::squish) +
  labs(x = "Astrocyte markers",
       y = expression(atop("Protein expression", "(" * Log[10] * " mean abundance)"))) +
  theme_pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
        axis.text.y = element_text(size = PT(SZ$axis + 0.5)),   # keeps the exponents legible
        legend.position = "none")
save_fig(p_prot, "Fig1A_astrocyte_markers", fig_w / 2, 2.6)


# ---- 5. Panel B: GFAP immunofluorescence ---------------------------------

read_image <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "png") png::readPNG(path)
  else if (ext %in% c("jpg","jpeg")) jpeg::readJPEG(path)
  else if (ext %in% c("tif","tiff")) {
    if (!requireNamespace("tiff", quietly = TRUE)) install.packages("tiff")
    tiff::readTIFF(path)
  } else stop("Unsupported image format: ", ext)
}

if (!is.na(IF_IMAGE) && file.exists(IF_IMAGE)) {
  img <- read_image(IF_IMAGE)
  p_if <- ggplot() +
    annotation_custom(rasterGrob(img, interpolate = TRUE,
                                 width = unit(1, "npc"), height = unit(1, "npc"))) +
    coord_fixed(ratio = dim(img)[1] / dim(img)[2], xlim = c(0, 1), ylim = c(0, 1),
                expand = FALSE) +
    theme_void()
  message("Panel B: using ", basename(IF_IMAGE))
} else {
  p_if <- ggplot() +
    annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1,
             fill = "grey97", colour = "grey60", linewidth = 0.4) +
    annotate("text", x = 0.5, y = 0.5, size = GEOM_TEXT, colour = "grey45",
             label = "GFAP immunofluorescence\n(set IF_IMAGE to the image file)") +
    coord_fixed(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) + theme_void()
  message("Panel B: image not found - placeholder used. Set IF_IMAGE.")
}
save_fig(p_if, "Fig1B_GFAP_immunofluorescence", fig_w / 2, 2.6)


# ---- 6. Read transcriptomic results --------------------------------------

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
  d$cls <- ifelse(!d$sig, "ns", ifelse(d$lfc > 0, "up", "down"))
  message(sprintf("  %-14s %5d DETs, %5d DEGs  <- %s", lab, sum(d$sig),
                  dplyr::n_distinct(d$gene[d$sig & !is.na(d$gene)]), basename(fp)))
  d
}

ALL <- list()
for (ln in names(LINES)) for (dg in names(DRUGS)) {
  if (any(vapply(EXCLUDE, function(e) e[1] == ln && e[2] == dg, logical(1)))) {
    message(sprintf("  %-14s excluded", paste(ln, dg))); next
  }
  d <- read_full(ln, dg); if (!is.null(d)) ALL[[paste(ln, dg)]] <- d
}
ALL <- bind_rows(ALL)
if (!nrow(ALL)) stop("No CSVs read - check analysis_dir and filename keywords.")
ALL <- ALL %>% mutate(drug = factor(drug, levels = names(DRUGS)),
                      line = factor(line, levels = names(LINES)))

SETS <- list()
for (ln in names(LINES)) for (dg in names(DRUGS)) {
  s <- ALL %>% filter(line == ln, drug == dg, sig)
  if (nrow(s)) SETS[[ln]][[dg]] <- s
}


# ---- 7. Panel C: volcano plots -------------------------------------------

x_lim <- c(floor(min(ALL$lfc)), ceiling(max(ALL$lfc)))
true_y_max <- max(-log10(ALL$padj), na.rm = TRUE)
y_disp <- if (is.na(Y_CAP)) ceiling(true_y_max) else Y_CAP

counts <- ALL %>% filter(sig) %>% group_by(drug, line) %>%
  summarise(n_up = sum(lfc > 0), n_dn = sum(lfc < 0), .groups = "drop")
ALLv <- ALL %>% mutate(y = pmin(-log10(padj), y_disp),
                       capped = -log10(padj) > y_disp,
                       cls = factor(cls, levels = c("ns","down","up")))

## cells with no data (SVG-A morphine) are labelled rather than left blank
blank <- expand_grid(line = levels(ALL$line), drug = levels(ALL$drug)) %>%
  anti_join(distinct(ALL, line, drug), by = c("line","drug")) %>%
  mutate(line = factor(line, levels = levels(ALL$line)),
         drug = factor(drug, levels = levels(ALL$drug)))

p_vol <- ggplot(ALLv, aes(lfc, y)) +
  geom_hline(yintercept = -log10(PADJ), linetype = "dashed", colour = "grey55", linewidth = 0.3) +
  annotate("segment", x = c(-LFC, LFC), xend = c(-LFC, LFC), y = 0, yend = y_disp * 1.0,
           linetype = "dashed", colour = "grey55", linewidth = 0.3) +
  geom_point(aes(colour = cls, shape = capped), size = 0.5, alpha = 0.6, stroke = 0.25) +
  ## counts in the top corners, coloured like the points they count
  geom_text(data = counts, aes(x = x_lim[1], y = y_disp * 1.12, label = paste(comma(n_dn), "down")),
            inherit.aes = FALSE, hjust = 0, vjust = 1, size = GEOM_TEXT,
            colour = COL_DOWN, fontface = "bold") +
  geom_text(data = counts, aes(x = x_lim[2], y = y_disp * 1.12, label = paste(comma(n_up), "up")),
            inherit.aes = FALSE, hjust = 1, vjust = 1, size = GEOM_TEXT,
            colour = COL_UP, fontface = "bold") +
  geom_text(data = blank, aes(x = 0, y = y_disp / 2, label = "excluded\n(flowcell confound)"),
            inherit.aes = FALSE, size = GEOM_TEXT, colour = "grey55", fontface = "italic") +
  scale_colour_manual(values = c(ns = COL_NS, down = COL_DOWN, up = COL_UP), guide = "none") +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 2), guide = "none") +
  scale_x_continuous(breaks = c(-20, -10, 0, 10, 20)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(xlim = x_lim, ylim = c(0, y_disp * 1.12)) +
  facet_grid(line ~ drug) +
  labs(x = expression(log[2]~"fold-change"), y = expression(-log[10]~adjusted~italic(p)),
       caption = sprintf("%s; grey, not significant. Shared axes; y capped at %g (open triangles); true max = %.0f.",
                         thr_lab, y_disp, true_y_max)) +
  theme_pub
save_fig(p_vol, "Fig1C_volcanoes", fig_w, 3.3)


# ---- 8. Panels D and E: cross-line Venn diagrams -------------------------

make_venn2 <- function(setK, setS, ttl) {
  nK <- length(setdiff(setK, setS)); nS <- length(setdiff(setS, setK))
  nB <- length(intersect(setK, setS))
  r <- 1.0; off <- 0.62; th <- seq(0, 2*pi, length.out = 200)
  circ <- function(cx, s) data.frame(x = cx + r*cos(th), y = r*sin(th), set = s)
  circles <- rbind(circ(-off, "KNS-42"), circ(off, "SVG-A"))
  ggplot() +
    geom_polygon(data = circles, aes(x, y, group = set, fill = set), alpha = 0.45, colour = NA) +
    geom_path(data = circles, aes(x, y, group = set, colour = set), linewidth = 0.4) +
    scale_fill_manual(values = setNames(VENN2_COLS, c("KNS-42","SVG-A")), guide = "none") +
    scale_colour_manual(values = setNames(VENN2_COLS, c("KNS-42","SVG-A")), guide = "none") +
    annotate("text", x = -off-0.45, y = 0, label = comma(nK), size = PT(SZ$venn_set) / .pt) +
    annotate("text", x =  off+0.45, y = 0, label = comma(nS), size = PT(SZ$venn_set) / .pt) +
    annotate("text", x = 0,         y = 0, label = comma(nB), size = PT(SZ$venn_set) / .pt) +
    annotate("text", x = -off, y = r + 0.08, label = "KNS-42", size = GEOM_TEXT, vjust = 0,
             fontface = "bold", colour = VENN2_COLS[1]) +
    annotate("text", x =  off, y = r + 0.08, label = "SVG-A", size = GEOM_TEXT, vjust = 0,
             fontface = "bold", colour = VENN2_COLS[2]) +
    labs(title = ttl) + coord_equal(clip = "off") +
    theme_void(base_size = PT(SZ$axis)) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = PT(SZ$strip),
                                    margin = margin(b = 12)),
          plot.margin = margin(4, 6, 4, 6))
}

venn_row <- function(level = c("transcript","gene")) {
  level <- match.arg(level)
  ps <- list(); summ <- list()
  for (dg in VENN_DRUGS) {
    k <- SETS[["KNS-42"]][[dg]]; s <- SETS[["SVG-A"]][[dg]]
    if (is.null(k) || is.null(s)) next
    if (level == "transcript") { sk <- k$id; ss <- s$id; ttl <- dg }
    else { sk <- unique(na.omit(k$gene)); ss <- unique(na.omit(s$gene)); ttl <- paste0(dg, " (genes)") }
    ps[[dg]] <- make_venn2(sk, ss, ttl)
    summ[[dg]] <- data.frame(drug = dg, level = level, KNS = length(sk), SVG = length(ss),
                             shared = length(intersect(sk, ss)))
  }
  list(plot = wrap_plots(ps, nrow = 1), summary = bind_rows(summ))
}

v_det <- venn_row("transcript"); v_deg <- venn_row("gene")
print(bind_rows(v_det$summary, v_deg$summary), row.names = FALSE)
write_csv(bind_rows(v_det$summary, v_deg$summary), file.path(out_dir, "Fig1_venn_counts.csv"))

p_det <- v_det$plot + plot_annotation(
  title = "DET overlap between cell lines, by drug (transcript level)", theme = title_theme)
p_deg <- v_deg$plot + plot_annotation(
  title = "DEG overlap between cell lines, by drug (gene level)", theme = title_theme)
save_fig(p_det, "Fig1D_DET_venn", fig_w, 1.55)
save_fig(p_deg, "Fig1E_DEG_venn", fig_w, 1.55)



# ---- 9. Assemble Figure 1 ------------------------------------------------

## A and B side by side on the top row, then C, D, E full width.
fig1 <- (wrap_elements(full = p_prot) | wrap_elements(full = p_if)) /
  wrap_elements(full = p_vol) /
  wrap_elements(full = p_det) /
  wrap_elements(full = p_deg) +
  plot_layout(heights = c(2.6, 3.3, 1.55, 1.55)) +
  plot_annotation(tag_levels = list(c("A", "B", "C", "D", "E"))) &
  tag_theme
save_fig(fig1, "Fig1_combined", fig_w, 9.0)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
message("Figure 1 written to: ", out_dir)
