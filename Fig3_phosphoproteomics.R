
# ---- 1. Settings ----------------------------------------------------------

analysis_dir <- "C:/Users/mrzsjs1/Desktop/Addiction/analysis"
out_dir      <- file.path(analysis_dir, "figures", "Fig3")

PHOS_DIR   <- "C:/Users/mrzsjs1/Downloads"          # *_phosphoproteomics_significant.xlsx
KINASE_DIR <- file.path(analysis_dir, "Kinase")      # *_Kinase_Scores.xlsx, *_Kinase_Enrichment_data.xlsx
SCORE_FILE  <- function(drug) file.path(KINASE_DIR, paste0(drug, "_Kinase_Scores.xlsx"))
ENRICH_FILE <- function(drug) file.path(KINASE_DIR, paste0(drug, "_Kinase_Enrichment_data.xlsx"))

KSEA_CSV <- file.path(KINASE_DIR, "kinase_figs", "kinase_activity_corrected_all.csv")
PTM_TSV  <- file.path("C:/Users/mrzsjs1/Downloads", "20241018_112329_P0904_Phospho_Report_PTM.tsv")

DRUGS <- c(Ethanol = "ethanol", THC = "thc", Cocaine = "cocaine",
           Morphine = "morphine", Psilocin = "psilocin")      # display order

PADJ <- 0.05; LFC <- 1
PHOS_P         <- 0.05     # phosphosites: raw p (few survive correction)
TOP_N_SITES    <- 10       # panel B
LABEL_BY_SITE  <- TRUE     # panel B labels "GENE S123" (FALSE = gene only)
VENN_NEEDS_LFC <- TRUE     # panel C: also require |log2FC| > 1 (as published)

## KSEA settings (recalculation only)
SITE_PROB <- 0.75; MIN_REP <- 3; MEDIAN_NORM <- TRUE
GROUPS   <- list(UC = 1:4, EC = 5:8, E = 9:12, T = 13:16, C = 17:20, M = 21:24, P = 25:28)
CONTRAST <- list(Ethanol = c("E","UC"), THC = c("T","EC"), Cocaine = c("C","UC"),
                 Morphine = c("M","UC"), Psilocin = c("P","UC"))
## heatmap settings
MIN_SUB <- 3; MIN_DRUGS <- 4; TOP_N <- 28; Z_CAP <- 4

# palette
COL_UP <- "#E24B4A"; COL_DOWN <- "#378ADD"
VENN_LOW <- "#F7FBFF"; VENN_HIGH <- "#7F77DD"
COL_INH <- "#2166AC"; COL_MID <- "#FFFFFF"; COL_ACT <- "#E08214"; NA_COL <- "grey96"


# ---- 2. Packages ----------------------------------------------------------

pkgs <- c("readr","readxl","dplyr","tidyr","tibble","ggplot2","scales","patchwork",
          "ggVennDiagram","data.table")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing)
suppressPackageStartupMessages({
  library(readr); library(readxl); library(dplyr); library(tidyr); library(ggplot2)
  library(scales); library(patchwork); library(ggVennDiagram); library(data.table)
})
select <- dplyr::select; filter <- dplyr::filter
rename <- dplyr::rename; mutate <- dplyr::mutate
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
num <- function(x) as.numeric(gsub(",", ".", x, fixed = TRUE))   # comma-decimal -> numeric

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

## kinase panels: kinases across, drugs down, legends underneath
theme_kin <- theme_minimal(base_size = PT(SZ$axis)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                   size = PT(SZ$annot), colour = "black"),
        axis.text.y = element_text(size = PT(SZ$tick), face = "bold", colour = "grey20"),
        panel.grid  = element_blank(),
        plot.title  = element_text(face = "bold", hjust = 0.5, size = PT(SZ$title),
                                   margin = margin(b = 3)),
        plot.title.position = "panel",
        legend.position = "bottom", legend.box = "horizontal",
        legend.title = element_text(size = PT(SZ$legend), vjust = 0.8),
        legend.text  = element_text(size = PT(SZ$legend)),
        legend.key.height = unit(6, "pt"), legend.key.width = unit(16, "pt"),
        legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(2, "pt"),
        plot.margin = margin(4, 6, 4, 4))
placeholder <- function(ttl, msg) {
  ggplot() + annotate("text", x = 0.5, y = 0.5, label = msg, size = GEOM_TEXT, colour = "grey45") +
    xlim(0, 1) + ylim(0, 1) + labs(title = ttl) + theme_void() + title_theme
}


# ---- 4. Read phosphoproteomics -------------------------------------------

find_phos <- function(key) {
  f <- list.files(PHOS_DIR, pattern = "\\.xlsx?$", full.names = TRUE)
  hit <- f[grepl(key, basename(f), ignore.case = TRUE) & grepl("phospho", basename(f), ignore.case = TRUE)]
  if (length(hit)) hit[1] else NA_character_
}
PHOS <- bind_rows(lapply(names(DRUGS), function(d) {
  fp <- find_phos(DRUGS[[d]])
  if (is.na(fp)) { warning("no phospho file for ", d); return(NULL) }
  read_excel(fp) %>%
    transmute(drug = d,
              gene = trimws(sub("[;,].*$", "", as.character(gene))),
              site = as.character(mod_site),
              lfc  = as.numeric(log2_fold_change),
              p    = as.numeric(p_value),
              padj = as.numeric(adj_p_value),
              ox   = grepl("Oxidation", as.character(sequence))) %>%
    filter(!is.na(gene), gene != "", !is.na(lfc), !is.na(p))   # unannotated sites dropped
})) %>% mutate(drug = factor(drug, levels = names(DRUGS)), sig = p < PHOS_P)

phos_counts <- PHOS %>% filter(sig) %>% group_by(drug) %>%
  summarise(n_up = sum(lfc > 0), n_dn = sum(lfc < 0),
            n_adj = sum(padj < PADJ, na.rm = TRUE), .groups = "drop")
print(as.data.frame(phos_counts), row.names = FALSE)
write_csv(phos_counts, file.path(out_dir, "Fig3A_phospho_counts.csv"))


# ---- 5. Panel A: volcano plots --------------------------------------------

PHv <- PHOS %>% filter(sig) %>%
  mutate(cls = factor(ifelse(lfc > 0, "Up", "Down"), levels = c("Up", "Down")))
xl   <- c(-1, 1) * ceiling(max(abs(PHv$lfc)))
ydat <- max(-log10(PHv$p))
yl   <- ceiling(ydat + 1.2)

p_A <- ggplot(PHv, aes(lfc, -log10(p))) +
  geom_hline(yintercept = -log10(PHOS_P), linetype = "dashed", colour = "grey55", linewidth = 0.3) +
  annotate("segment", x = c(-LFC, LFC), xend = c(-LFC, LFC), y = 0, yend = ydat + 0.4,
           linetype = "dashed", colour = "grey55", linewidth = 0.3) +
  geom_point(aes(colour = cls), size = 0.7, alpha = 0.7, stroke = 0) +
  geom_text(data = phos_counts, aes(x = xl[1], y = yl, label = paste(n_dn, "down")),
            inherit.aes = FALSE, hjust = 0, vjust = 1, size = GEOM_TEXT,
            colour = COL_DOWN, fontface = "bold") +
  geom_text(data = phos_counts, aes(x = xl[2], y = yl, label = paste(n_up, "up")),
            inherit.aes = FALSE, hjust = 1, vjust = 1, size = GEOM_TEXT,
            colour = COL_UP, fontface = "bold") +
  scale_colour_manual(values = c(Up = COL_UP, Down = COL_DOWN), guide = "none") +
  scale_x_continuous(breaks = c(-5, 0, 5)) +
  scale_y_continuous(breaks = seq(0, 20, 2), expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(xlim = xl, ylim = c(0, yl)) +
  facet_wrap(~ drug, nrow = 1) +
  labs(title = "Phosphoproteomic response by drug (KNS-42)",
       x = expression(log[2]~"fold-change"), y = expression(-log[10]~italic(p)),
       caption = paste0("Raw p < 0.05. Sites also significant after multiple-testing correction (adj. p < 0.05): ",
                        paste(sprintf("%s %d", phos_counts$drug, phos_counts$n_adj), collapse = ", "), ".")) +
  theme_pub
save_fig(p_A, "Fig3A_phospho_volcanoes", fig_w, 1.7)


# ---- 6. Panel B: top regulated phosphosites -------------------------------

top <- PHOS %>% filter(sig) %>%
  mutate(site = gsub(",\\s*", "/", site),
         label = if (LABEL_BY_SITE) paste(gene, site) else gene) %>%
  group_by(drug) %>%
  slice_max(abs(lfc), n = TOP_N_SITES, with_ties = FALSE) %>%
  ## the same site can be quantified on more than one peptide (e.g. THC TJP1 Y895,
  ## with and without Met oxidation): tag the oxidised form, number any others
  mutate(dup = duplicated(label) | duplicated(label, fromLast = TRUE),
         label = ifelse(dup & ox, paste(label, "(M-ox)"), label),
         label = ifelse(duplicated(label),
                        paste0(label, " (", ave(label, label, FUN = seq_along), ")"), label)) %>%
  arrange(lfc, .by_group = TRUE) %>%
  mutate(yk = factor(paste(label, drug, sep = "___"),
                     levels = unique(paste(label, drug, sep = "___")))) %>%
  ungroup() %>%
  mutate(cls = factor(ifelse(lfc > 0, "Up", "Down"), levels = c("Up", "Down")))
write_csv(top %>% select(drug, label, gene, site, lfc, p, padj, met_oxidised = ox),
          file.path(out_dir, "Fig3B_top_sites.csv"))

bx  <- ceiling(max(abs(top$lfc)))
drs <- levels(droplevels(top$drug))
one_lolli <- function(d, xlab = FALSE) {
  ggplot(dplyr::filter(top, drug == d), aes(lfc, yk, colour = cls)) +
    geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.4) +
    geom_segment(aes(x = 0, xend = lfc, yend = yk), linewidth = 0.4) +
    geom_point(size = 1.1) +
    scale_colour_manual(values = c(Up = COL_UP, Down = COL_DOWN), guide = "none") +
    scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
    scale_x_continuous(limits = c(-bx, bx), breaks = pretty_breaks(3)) +
    facet_wrap(~ drug) +
    labs(x = if (xlab) expression(log[2]~"fold-change") else NULL, y = NULL) +
    theme_pub +
    theme(axis.text.y = element_text(size = PT(SZ$annot)))
}
## two rows, the shorter bottom row centred under the top row
n1  <- ceiling(length(drs) / 2); n2 <- length(drs) - n1
row1 <- wrap_plots(lapply(drs[seq_len(n1)], one_lolli), nrow = 1)
pad  <- (n1 - n2) / 2
row2 <- wrap_plots(c(list(plot_spacer()),
                     lapply(drs[n1 + seq_len(n2)], one_lolli, xlab = TRUE),
                     list(plot_spacer())),
                   nrow = 1, widths = c(pad, rep(1, n2), pad))
p_B <- (row1 / row2) +
  plot_annotation(title = sprintf("Top %d regulated phosphosites per drug (KNS-42)", TOP_N_SITES),
                  theme = title_theme)
save_fig(p_B, "Fig3B_top_phosphosites", fig_w, 3.7)


# ---- 7. Panel C: phosphoproteins shared across drugs ----------------------

## As published, the Venn uses raw p < 0.05 AND |log2FC| > 1; panels A/B show
## every raw p < 0.05 site. VENN_NEEDS_LFC <- FALSE uses the A/B rule.
vsig <- PHOS$sig & (!VENN_NEEDS_LFC | abs(PHOS$lfc) > LFC)
phos_sets <- lapply(split(PHOS$gene[vsig], PHOS$drug[vsig]), unique)
phos_sets <- phos_sets[lengths(phos_sets) > 0]
p_C <- make_venn(phos_sets, "Shared phosphoproteins", "proteins") +
  theme(legend.position = "bottom", legend.title = element_text(size = PT(SZ$legend), vjust = 0.8),
        legend.key.height = unit(6, "pt"), legend.key.width = unit(16, "pt"))
venn_counts <- data.frame(drug = names(phos_sets), total = lengths(phos_sets),
  private = vapply(names(phos_sets), function(d)
    length(setdiff(phos_sets[[d]], unlist(phos_sets[names(phos_sets) != d]))), integer(1)),
  in_all = length(Reduce(intersect, phos_sets)))
print(venn_counts, row.names = FALSE)
write_csv(venn_counts, file.path(out_dir, "Fig3C_venn_counts.csv"))
save_fig(p_C, "Fig3C_phosphoprotein_venn", 2.8, 4.0)


# ---- 8. Panel D: kinase substrate over-representation ---------------------

KDRUGS <- names(CONTRAST)
miss <- KDRUGS[!file.exists(vapply(KDRUGS, ENRICH_FILE, character(1)))]
if (length(miss)) {
  message("Panel D skipped; missing enrichment files for: ", paste(miss, collapse = ", "))
  p_D <- placeholder("Kinase substrate enrichment", "Enrichment files not found\n(see KINASE_DIR)")
} else {
  enr <- bind_rows(lapply(KDRUGS, function(drug) {
    e <- read_excel(ENRICH_FILE(drug))
    tibble(kinase = as.character(e$Description), drug = drug,
           neglog10 = -log10(num(e$`p.adjust`)), count = num(e$Count))
  }))
  write_csv(enr, file.path(out_dir, "Fig3D_kinase_enrichment_all.csv"))
  topk <- enr %>% group_by(kinase) %>%
    summarise(s = max(neglog10, na.rm = TRUE), .groups = "drop") %>%
    slice_max(s, n = TOP_N) %>% arrange(desc(s)) %>% pull(kinase)  # most significant on the left
  p_D <- enr %>% filter(kinase %in% topk) %>%
    mutate(drug = factor(drug, levels = rev(KDRUGS)), kinase = factor(kinase, levels = topk)) %>%
    ggplot(aes(kinase, drug, size = count, colour = neglog10)) +
    geom_point() +
    scale_colour_gradient(low = "#F6D9BE", high = COL_ACT, name = "-log10 p.adj") +
    scale_size_area(max_size = 3.2, name = "Substrates", breaks = c(10, 20, 40)) +
    guides(colour = guide_colourbar(order = 1),
           size = guide_legend(order = 2, override.aes = list(colour = "black"))) +
    labs(title = "Kinase substrate enrichment", x = NULL, y = NULL) +
    theme_kin +
    theme(panel.grid.major = element_line(colour = "grey93", linewidth = 0.25))
}
save_fig(p_D, "Fig3D_kinase_enrichment", 4.3, 2.0)


# ---- 9. Panel E: predicted kinase activity (KSEA) -------------------------

use_csv <- file.exists(KSEA_CSV)
ksea_inputs <- c(PTM_TSV, vapply(KDRUGS, SCORE_FILE, character(1)))
if (!use_csv && !all(file.exists(ksea_inputs))) {
  message("Panel E skipped; set KSEA_CSV, or PTM_TSV and KINASE_DIR")
  p_E <- placeholder("Kinase activity (KSEA)", "KSEA inputs not found\n(set KSEA_CSV / PTM_TSV)")
} else {
  if (use_csv) {
    message("Panel E: saved full-background KSEA scores: ", KSEA_CSV)
    score <- read_csv(KSEA_CSV, show_col_types = FALSE) %>%
      transmute(kinase = as.character(kinase), z = as.numeric(z), p = as.numeric(p),
                m = as.numeric(m), drug = as.character(drug), FDR = as.numeric(FDR))
  } else {
    message("Reading PTM report ...")
    dt <- fread(PTM_TSV, sep = "\t",
                select = c("R.FileName","PTM.CollapseKey","PTM.ModificationTitle",
                           "PTM.Quantity","PTM.SiteProbability"),
                colClasses = "character", showProgress = FALSE)
    dt <- dt[PTM.ModificationTitle == "Phospho (STY)"]
    dt[, samp := as.integer(sub(".*_phospho_([0-9]+).*", "\\1", R.FileName))]
    dt[, q    := num(PTM.Quantity)]
    dt[, prob := num(PTM.SiteProbability)]
    dt <- dt[grepl("_[STY][0-9]+_M[0-9]+$", PTM.CollapseKey)]
    dt[, key := sub("_M[0-9]+$", "", PTM.CollapseKey)]
    dt <- dt[!is.na(q) & q > 0 & !is.na(prob) & prob >= SITE_PROB & !is.na(samp)]
    agg  <- dt[, .(q = sum(q)), by = .(key, samp)]
    wide <- dcast(agg, key ~ samp, value.var = "q", fill = NA)
    M <- as.matrix(log2(wide[, -1])); rownames(M) <- wide$key
    colnames(M) <- as.integer(colnames(M))
    if (MEDIAN_NORM) M <- sweep(M, 2, apply(M, 2, median, na.rm = TRUE), "-")
    bg <- do.call(rbind, lapply(CONTRAST, function(cc) {
      dcol <- as.character(GROUPS[[cc[1]]]); ccol <- as.character(GROUPS[[cc[2]]])
      keep <- rowSums(!is.na(M[, dcol])) >= MIN_REP & rowSums(!is.na(M[, ccol])) >= MIN_REP
      fc <- rowMeans(M[keep, dcol], na.rm = TRUE) - rowMeans(M[keep, ccol], na.rm = TRUE)
      c(mP = mean(fc), sdP = sd(fc), n = sum(keep))
    }))
    rownames(bg) <- KDRUGS
    print(round(bg, 3))
    score <- bind_rows(lapply(KDRUGS, function(drug) {
      s <- read_excel(SCORE_FILE(drug)); s$mS <- num(s$mS); s$m <- num(s$m)
      mP <- bg[drug, "mP"]; sdP <- bg[drug, "sdP"]
      s %>% filter(!is.na(mS), !is.na(m)) %>%
        transmute(kinase = Kinase.Gene, z = (mS - mP) * sqrt(m) / sdP,
                  p = 2 * pnorm(-abs(z)), m = m, drug = drug)
    })) %>% group_by(drug) %>% mutate(FDR = p.adjust(p, "BH")) %>% ungroup()
  }
  write_csv(score, file.path(out_dir, "Fig3E_kinase_activity_all.csv"))

  wide_of <- function(col) {
    w <- score %>% filter(m >= MIN_SUB) %>%
      pivot_wider(id_cols = kinase, names_from = drug, values_from = all_of(col)) %>%
      as.data.frame()
    rownames(w) <- w$kinase; w$kinase <- NULL; w[, KDRUGS]
  }
  zw <- wide_of("z"); fdrw <- wide_of("FDR")
  zw <- zw[rowSums(!is.na(zw)) >= MIN_DRUGS, , drop = FALSE]
  fdrw <- fdrw[rownames(zw), , drop = FALSE]

  ## kinases with identical z profiles (family aliases) collapsed into one row
  keyvec <- apply(round(zw, 3), 1, function(r) paste(ifelse(is.na(r), "NA", r), collapse = "|"))
  reps  <- tapply(rownames(zw), keyvec, function(v) {
    v <- sort(v); paste0(paste(head(v, 3), collapse = "/"), if (length(v) > 3) "\u2026" else "")
  })
  first <- tapply(rownames(zw), keyvec, function(v) sort(v)[1])
  zc <- zw[first, , drop = FALSE];   rownames(zc) <- reps[names(first)]
  fc <- fdrw[first, , drop = FALSE]; rownames(fc) <- reps[names(first)]
  ord <- order(rowMeans(abs(zc), na.rm = TRUE), decreasing = TRUE)
  zc <- zc[head(ord, TOP_N), , drop = FALSE]
  zc <- zc[order(rowMeans(zc, na.rm = TRUE), decreasing = TRUE), , drop = FALSE]  # most activated on the left
  fc <- fc[rownames(zc), , drop = FALSE]

  klong <- zc %>% tibble::rownames_to_column("kinase") %>%
    pivot_longer(-kinase, names_to = "drug", values_to = "z") %>%
    left_join(fc %>% tibble::rownames_to_column("kinase") %>%
                pivot_longer(-kinase, names_to = "drug", values_to = "FDR"),
              by = c("kinase","drug")) %>%
    mutate(drug = factor(drug, levels = rev(KDRUGS)),
           kinase = factor(kinase, levels = rownames(zc)),
           z = pmax(pmin(z, Z_CAP), -Z_CAP))
  write_csv(klong, file.path(out_dir, "Fig3E_kinase_activity_shown.csv"))

  p_E <- ggplot(klong, aes(kinase, drug, fill = z)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_point(data = subset(klong, !is.na(FDR) & FDR < 0.05),
               shape = 21, fill = "white", colour = "grey15", size = 0.9, stroke = 0.3) +
    scale_fill_gradient2(low = COL_INH, mid = COL_MID, high = COL_ACT, midpoint = 0,
                         na.value = NA_COL, limits = c(-Z_CAP, Z_CAP), name = "KSEA z-score") +
    labs(title = "Kinase activity (KSEA)", x = NULL, y = NULL) +
    theme_kin
}
save_fig(p_E, "Fig3E_kinase_activity_KSEA", 4.3, 2.0)


# ---- 10. Assemble Figure 3 -------------------------------------------------

bottom <- (wrap_elements(full = p_C) | (p_D / p_E)) + plot_layout(widths = c(2.8, 4.3))
fig3 <- wrap_elements(full = p_A) /
  wrap_elements(full = p_B) /
  bottom +
  plot_layout(heights = c(1.7, 3.7, 4.0)) +
  plot_annotation(tag_levels = list(c("A", "B", "C", "D", "E"))) &
  tag_theme
save_fig(fig3, "Fig3_combined", fig_w, 9.4)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
message("Figure 3 written to: ", out_dir)
