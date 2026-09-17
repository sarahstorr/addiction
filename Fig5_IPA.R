
pkgs <- c("data.table","dplyr","tidyr","stringr","ggplot2","scales","patchwork")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressPackageStartupMessages(invisible(lapply(pkgs, library, character.only = TRUE)))
select <- dplyr::select; filter <- dplyr::filter; slice <- dplyr::slice

## ============================ FILES ======================================
IPA_DIR   <- "C:/Users/mrzsjs1/Desktop/Addiction/reassessIPA"
OUT_DIR   <- file.path(IPA_DIR, "Fig5")
out_dir   <- OUT_DIR
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

F_CANON_KNS   <- "txnphos_canonical.txt"   # panel B
F_UPSTR_KNS   <- "txnphos_upstream.txt"    # panel D
F_CANON_CROSS <- "lines_cannonical.txt"    # panel A
F_UPSTR_CROSS <- "lines_upstream.txt"      # panel C

DRUGS_ALL    <- c("Ethanol","THC","Cocaine","Morphine","Psilocin")   # within-KNS-42
DRUGS_NOMORP <- c("Ethanol","THC","Cocaine","Psilocin")              # cross-line (no morphine)

# ---- Selection rules (applied identically to all four panels) -------------------

## Ranking: how rows are ordered and how many are kept
SELECT_MODE <- "recurrence"  # "recurrence" = rank by the number of comparisons
                             #   in which a z-score could be computed, ties
                             #   broken by mean absolute z
                             # "magnitude"  = rank by mean absolute z only
TOP_N_PATH  <- 25            # canonical pathways displayed
TOP_N_REG   <- 24            # upstream regulators displayed

ROW_ORDER   <- "cluster"     # how the SELECTED rows are ordered top-to-bottom.
                             # "cluster" = hierarchical clustering (average
                             #   linkage on 1 - pairwise-complete correlation),
                             #   so co-regulated pathways sit together
                             # "rank"    = ranked order (recurrence/magnitude)
                             # (affects order only, not which rows are shown)

## Eligibility: which rows may enter the ranking
MIN_LAYERS    <- 2      # a row must score in at least this many of the two
                        #   layers compared (2 = scored in both layers)
DROP_ALL_ZERO <- TRUE   # drop rows whose z-scores are all exactly 0

## Row types
EXCLUDE_SUPERPATHWAYS <- TRUE   # pathway panels: drop broad disease terms
EXCLUDE_PATH <- c("Cardiac Hypertrophy","Hepatic Fibrosis","Molecular Mechanisms of Cancer",
                  "Colorectal Cancer","Chronic Myeloid Leukemia","Autism","Atherosclerosis",
                  "Rheumatoid","Osteoarthritis","Tumor Microenvironment","Pancreatic",
                  "Glioblastoma","Breast Cancer","Systemic Lupus","Pulmonary Fibrosis")

GENELIKE_ONLY    <- TRUE   # regulator panels: gene/protein symbols only
                           #   (drops chemicals, drugs and microRNAs)
GENELIKE_MAXCHAR <- 8      # max symbol length kept when GENELIKE_ONLY is TRUE.
                           #   (raise to ~15 to keep long/fusion symbols)

## Compounds whose names resemble gene symbols (not removed by GENELIKE_ONLY)
DROP_CHEMICALS <- TRUE
CHEMICAL_NAMES <- c("PD98059",   # MEK1/2 inhibitor
                    "SP2509",    # LSD1/KDM1A inhibitor
                    "ST1926",    # atypical retinoid
                    "U0126", "LY294002", "SB203580", "SP600125", "PD0325901",
                    "GW4064", "T0901317", "AICAR", "TSA", "SAHA")

## Display
Z_CAP   <- 4               # z-scores clamped to +/- this for the colour scale
                           #   and for the mean |z| tie-break
COL_INH <- "#2166AC"; COL_MID <- "#FFFFFF"; COL_ACT <- "#E08214"; NA_COL <- "grey93"

## Optional fixed row lists (override the ranking; names as in column 1 of the IPA file)
PIN_CANON_KNS   <- NULL   # panel B
PIN_UPSTR_KNS   <- NULL   # panel D
PIN_CANON_CROSS <- NULL   # panel A
PIN_UPSTR_CROSS <- NULL   # panel C

STRICT_COLUMNS <- TRUE     # stop if any IPA column cannot be mapped to a layer and drug

print_selection_rules <- function(){
  message("\n", strrep("=", 74))
  message("Selection rules in effect (all four panels)")
  message(strrep("=", 74))
  message("  ranking            : ", if (SELECT_MODE == "recurrence")
          "recurrence (n comparisons with a z-score), ties by mean |z|" else
          "mean |z| only")
  message("  displayed          : top ", TOP_N_PATH, " pathways / ", TOP_N_REG, " regulators")
  message("  row order          : ", if (ROW_ORDER == "cluster")
          "hierarchical clustering (average linkage, pairwise-complete corr)" else
          "ranked order")
  message("  min layers scored  : ", MIN_LAYERS, " of 2",
          if (MIN_LAYERS >= 2) "  (row must appear in BOTH panels)" else "  (no cross-layer requirement)")
  message("  all-zero rows      : ", if (DROP_ALL_ZERO) "dropped" else "kept")
  message("  super-pathways     : ", if (EXCLUDE_SUPERPATHWAYS)
          paste0("excluded (", length(EXCLUDE_PATH), " terms)") else "kept")
  message("  regulator rows     : ", if (GENELIKE_ONLY)
          paste0("gene symbols only, <= ", GENELIKE_MAXCHAR,
                 " chars (chemicals/drugs/miRNAs excluded)") else "all regulator types kept")
  message("  z-score cap        : +/- ", Z_CAP)
  message(strrep("=", 74), "\n")
}

CANON <- c(ethanol="Ethanol", thc="THC", cocaine="Cocaine", morphine="Morphine", psilocin="Psilocin")
LAYER_LABELS <- c(TXN_KNS="KNS-42", TXN_SVG="SVG-A", PHOS="Phosphoproteomics", METAB="Metabolomics")
LAYER_ORDER  <- c("TXN_KNS","TXN_SVG","PHOS","METAB")

## ==================== readers / parsers ==================================
read_ipa <- function(path){
  if (!file.exists(path)) stop("File not found: ", path)
  raw <- readLines(path, warn = FALSE, encoding = "latin1")
  hdr <- which(grepl("^(Upstream Regulators|Canonical Pathways|Ingenuity Canonical Pathways)\t", raw))[1]
  if (is.na(hdr)) hdr <- which(grepl("Upstream Regulators|Canonical Pathways", raw))[1]
  if (is.na(hdr)) stop("No IPA header row found in ", basename(path))
  d <- fread(path, sep = "\t", skip = hdr - 1, header = TRUE, encoding = "Latin-1", na.strings = "N/A")
  setnames(d, 1, "reg")
  setnames(d, names(d)[-1], sub("\\s*-\\s*20\\d\\d.*$", "", names(d)[-1]))
  d[, reg := trimws(as.character(reg))]
  d
}
layer_of <- function(x)
  ifelse(grepl("^IPA_metab",   x), "METAB",
  ifelse(grepl("^IPA_phospho", x), "PHOS",
  ifelse(grepl("^IPA_KNS",     x), "TXN_KNS",
  ifelse(grepl("^IPA_SVG",     x), "TXN_SVG", NA))))

drug_of <- function(x){
  s <- sub("^IPA_(metab_HMDB|phospho|KNS|SVG)_", "", x)
  s <- sub("\\s*\\(.*\\)\\s*$", "", s)      # drop "(reassess)" and similar tags
  s <- sub("[_ -]+(v\\d+|reassess|corrected|new|final)$", "", s, ignore.case = TRUE)
  unname(CANON[tolower(trimws(s))])
}

num <- function(v) suppressWarnings(as.numeric(trimws(gsub("[\r\"]", "", as.character(v)))))

to_long <- function(d){
  cols <- names(d)[-1]
  bind_rows(lapply(cols, function(c)
    tibble(reg = d$reg, layer = layer_of(c), drug = drug_of(c), z = num(d[[c]]), column = c)))
}

report_ipa <- function(long, name){
  message("-- ", name)
  map <- long %>% distinct(column, layer, drug) %>%
    left_join(long %>% group_by(column) %>% summarise(n = sum(!is.na(z)), .groups = "drop"), by = "column")
  for (i in seq_len(nrow(map))) {
    bad <- is.na(map$layer[i]) || is.na(map$drug[i])
    message(sprintf("   %-34s %-8s %-9s n=%-5d%s", map$column[i],
                    ifelse(is.na(map$layer[i]), "??", map$layer[i]),
                    ifelse(is.na(map$drug[i]),  "??", map$drug[i]),
                    map$n[i], if (bad) "   <- unmapped, would be dropped" else ""))
  }
  bad <- map %>% filter(is.na(layer) | is.na(drug))
  if (nrow(bad)) {
    msg <- paste0(nrow(bad), " column(s) in ", name, " could not be mapped: ",
                  paste(bad$column, collapse = " | "))
    if (STRICT_COLUMNS) stop(msg, "\n  -> fix layer_of()/drug_of() or rename the IPA analysis.")
    else warning(msg)
  }
  invisible(long)
}

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
select_ipa <- function(long, outstem,
                     kind = c("regulator","pathway"),
                     drugs = DRUGS_ALL, layers = "auto", layer_labels = NULL,
                     pin = NULL){
  kind <- match.arg(kind)
  if (identical(layers, "auto")) layers <- LAYER_ORDER[LAYER_ORDER %in% unique(long$layer[!is.na(long$layer)])]
  layers <- layers[layers %in% unique(long$layer[!is.na(long$layer)])]
  if (!length(layers)) stop("[", basename(outstem), "] requested layers not present.")
  if (is.null(layer_labels)) layer_labels <- LAYER_LABELS
  top_n <- if (kind == "pathway") TOP_N_PATH else TOP_N_REG

  L <- long %>% filter(layer %in% layers, drug %in% drugs)

  ## row type filter
  if (kind == "regulator" && GENELIKE_ONLY) {
    L <- L %>% filter(reg == toupper(reg), grepl("[A-Z]", reg),
                      !grepl(" ", reg), nchar(reg) <= GENELIKE_MAXCHAR)
    if (DROP_CHEMICALS) {
      hit <- intersect(unique(L$reg), CHEMICAL_NAMES)
      if (length(hit))
        message("   excluded chemical regulators: ", paste(hit, collapse = ", "))
      L <- L %>% filter(!reg %in% CHEMICAL_NAMES)
    }
  }
  if (!nrow(L %>% filter(!is.na(z))))
    stop("[", basename(outstem), "] no non-NA z after filtering. See report line above.")

  if (!is.null(pin)) {
    present <- pin[pin %in% L$reg]
    missing <- setdiff(pin, present)
    if (length(missing))
      warning("[", basename(outstem), "] ", length(missing),
              " pinned row(s) not found and skipped: ", paste(missing, collapse = " | "))
    if (!length(present)) stop("[", basename(outstem), "] none of the pinned rows were found.")
    sel    <- tibble(reg = rev(present))        # first-listed appears at top
    n_elig <- length(present)
  } else {
    Lf <- L %>% filter(!is.na(z))
    if (kind == "pathway" && EXCLUDE_SUPERPATHWAYS && length(EXCLUDE_PATH))
      Lf <- Lf %>% filter(!grepl(paste(EXCLUDE_PATH, collapse = "|"), reg, ignore.case = TRUE))

    stat <- Lf %>% group_by(reg) %>%
      summarise(nl  = n_distinct(layer),
                cov = n(),                                            # recurrence
                s   = mean(abs(pmin(pmax(z, -Z_CAP), Z_CAP))),        # tie-break
                .groups = "drop")

    ## eligibility
    stat <- stat %>% filter(nl >= min(MIN_LAYERS, length(layers)))
    if (DROP_ALL_ZERO) stat <- stat %>% filter(s > 0)
    n_elig <- nrow(stat)

    stat <- if (SELECT_MODE == "recurrence") stat %>% arrange(desc(cov), desc(s)) else
                                             stat %>% arrange(desc(s))
    sel <- stat %>% slice_head(n = top_n) %>% arrange(if (SELECT_MODE == "recurrence") cov else s)
    if (!nrow(sel)) stop("[", basename(outstem), "] nothing passed selection.")
  }

  ## optional hierarchical clustering of the selected rows
  ## Ordering only: the row set is already fixed by the selection rules. Distances use
  ## pairwise-complete correlation, so blank cells are skipped rather than
  ## treated as zero (an unscored pathway is "not assessed", not "no change").
  if (ROW_ORDER == "cluster" && nrow(sel) > 2) {
    M <- L %>% filter(reg %in% sel$reg) %>%
      mutate(cell = paste(layer, drug, sep = "_")) %>%
      select(reg, cell, z) %>%
      pivot_wider(names_from = cell, values_from = z,
                  values_fn = function(x) mean(x, na.rm = TRUE))
    Mm <- as.matrix(M[, -1, drop = FALSE]); rownames(Mm) <- M$reg
    keep <- sel$reg[sel$reg %in% rownames(Mm)]
    Mm <- Mm[keep, , drop = FALSE]
    if (nrow(Mm) > 2) {
      C <- suppressWarnings(stats::cor(t(Mm), use = "pairwise.complete.obs"))
      C[!is.finite(C)] <- 0
      D <- 1 - C; D[!is.finite(D)] <- 1; diag(D) <- 0
      hc  <- stats::hclust(stats::as.dist(D), method = "average")
      ord <- rownames(Mm)[hc$order]
      sel <- sel %>% slice(match(ord, sel$reg))
      message("   rows ordered by hierarchical clustering")
    }
  }

  ylab <- sel$reg   # single-line labels (print-size layout)
  P <- expand_grid(reg = sel$reg, drug = drugs, layer = layers) %>%
    left_join(L %>% select(reg, drug, layer, z), by = c("reg","drug","layer")) %>%
    mutate(reg   = factor(reg, levels = sel$reg, labels = ylab),
           drug  = factor(drug, levels = drugs),
           layer = factor(layer, levels = layers, labels = layer_labels[layers]),
           z     = pmax(pmin(z, Z_CAP), -Z_CAP))

  write.csv(P, paste0(outstem, "_matrix.csv"), row.names = FALSE)

  message("   ", nrow(sel), " of ", n_elig, " eligible ", kind, "s shown -> ", basename(outstem))
  P
}

# ---- Heatmap panel ------------------------------------------------------------

theme_heat <- theme_minimal(base_size = PT(SZ$axis)) +
  theme(axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                    size = PT(SZ$annot), colour = "grey20"),
        axis.text.y  = element_text(size = PT(SZ$annot), colour = "black"),
        panel.grid   = element_blank(),
        strip.text   = element_text(face = "bold", size = PT(SZ$annot), margin = margin(2, 1, 2, 1)),
        plot.title   = element_text(face = "bold", hjust = 0.5, size = PT(SZ$title),
                                    margin = margin(b = 3)),
        plot.title.position = "panel",
        panel.spacing = unit(9, "pt"),
        strip.clip = "off",                 # long layer names may extend past the strip
        legend.position = "bottom",
        legend.title = element_text(size = PT(SZ$legend), vjust = 0.8),
        legend.text  = element_text(size = PT(SZ$legend)),
        legend.key.height = unit(6, "pt"), legend.key.width = unit(22, "pt"),
        plot.margin = margin(4, 6, 4, 4))

heat_panel <- function(P, title, right_margin = 6) {
  ggplot(P, aes(drug, reg, fill = z)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    scale_fill_gradient2(low = COL_INH, mid = COL_MID, high = COL_ACT, midpoint = 0,
                         na.value = NA_COL, limits = c(-Z_CAP, Z_CAP),
                         breaks = seq(-Z_CAP, Z_CAP, 2), name = "Activation z-score") +
    scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
    facet_grid(~ layer) +
    labs(title = title, x = NULL, y = NULL) +
    theme_heat +
    theme(plot.margin = margin(4, right_margin, 4, 4))
}


# ---- Panels ---------------------------------------------------------------------

print_selection_rules()

## A. Canonical pathways - cross-line, morphine excluded
cp_nm <- to_long(read_ipa(file.path(IPA_DIR, F_CANON_CROSS)))
report_ipa(cp_nm, "A. Canonical pathways, cross-line (no morphine)")
P_A <- select_ipa(cp_nm, file.path(OUT_DIR, "Fig5A_canonical_crossline"),
                  kind = "pathway", drugs = DRUGS_NOMORP,
                  layers = c("TXN_KNS","TXN_SVG"),
                  layer_labels = c(TXN_KNS = "KNS-42", TXN_SVG = "SVG-A"),
                  pin = PIN_CANON_CROSS)

## B. Canonical pathways - KNS-42 phospho + transcriptomics
cp_kns <- to_long(read_ipa(file.path(IPA_DIR, F_CANON_KNS)))
report_ipa(cp_kns, "B. Canonical pathways, KNS-42 phospho + txn")
P_B <- select_ipa(cp_kns, file.path(OUT_DIR, "Fig5B_canonical_kns"),
                  kind = "pathway", drugs = DRUGS_ALL,
                  layers = c("PHOS","TXN_KNS"),
                  layer_labels = c(PHOS = "Phosphoproteome", TXN_KNS = "Transcriptome"),
                  pin = PIN_CANON_KNS)

## C. Upstream regulators - cross-line, morphine excluded
ur_nm <- to_long(read_ipa(file.path(IPA_DIR, F_UPSTR_CROSS)))
report_ipa(ur_nm, "C. Upstream regulators, cross-line (no morphine)")
P_C <- select_ipa(ur_nm, file.path(OUT_DIR, "Fig5C_upstream_crossline"),
                  kind = "regulator", drugs = DRUGS_NOMORP,
                  layers = c("TXN_KNS","TXN_SVG"),
                  layer_labels = c(TXN_KNS = "KNS-42", TXN_SVG = "SVG-A"),
                  pin = PIN_UPSTR_CROSS)

## D. Upstream regulators - KNS-42 phospho + transcriptomics
ur_kns <- to_long(read_ipa(file.path(IPA_DIR, F_UPSTR_KNS)))
report_ipa(ur_kns, "D. Upstream regulators, KNS-42 phospho + txn")
P_D <- select_ipa(ur_kns, file.path(OUT_DIR, "Fig5D_upstream_kns"),
                  kind = "regulator", drugs = DRUGS_ALL,
                  layers = c("PHOS","TXN_KNS"),
                  layer_labels = c(PHOS = "Phosphoproteome", TXN_KNS = "Transcriptome"),
                  pin = PIN_UPSTR_KNS)

p_A <- heat_panel(P_A, "Canonical pathways: cell lines")
p_B <- heat_panel(P_B, "Canonical pathways: KNS-42")
p_C <- heat_panel(P_C, "Upstream regulators: cell lines", right_margin = 10)
p_D <- heat_panel(P_D, "Upstream regulators: KNS-42",     right_margin = 10)

save_fig(p_A, "Fig5A_canonical_crossline", 4.3, 4.2)
save_fig(p_B, "Fig5B_canonical_kns",       4.3, 4.2)
save_fig(p_C, "Fig5C_upstream_crossline",  2.8, 4.1)
save_fig(p_D, "Fig5D_upstream_kns",        2.8, 4.1)


# ---- Assemble Figure 5 ------------------------------------------------------------
# Pathway panels (long names) stacked on the left, regulator panels on the right.

fig5 <- ((p_A / p_B) | (p_C / p_D)) +
  ## widths refer to the heatmap areas (row names are added on top)
  plot_layout(widths = c(2.1, 2.0), guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom") & tag_theme
save_fig(fig5, "Fig5_combined", fig_w, 8.6)

message("\nDONE. Figure 5 in: ", OUT_DIR)
