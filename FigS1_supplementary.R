## ========================================================================= ##
##  Supplementary Figure S1 - baseline target expression, transcriptional
##  assessment of cytotoxicity and additional drug-specific analyses
##
##  (A) Baseline expression of canonical drug targets (media-only controls)
##  (B) Cytotoxicity and stress hallmark genes at 24 h
##  (C) Drug-unique genes across all drugs
##  (D) Psilocin-unique genes in both cell lines
##  (E) Psilocin-unique metabolite features
##
##  Self-contained: every panel is computed here. Also writes
##    TableS1_baseline_targets.csv      baseline expression, all target genes
##    TableS2_cytotoxicity_sets.csv     hallmark set tests
##    TableS2_cytotoxicity_genes.csv    hallmark genes, per treatment
##    baseline_gene_TPM.csv / _counts   all genes, control samples
##    baseline_text_summary.txt         one line per target gene
##    FigS1_legend.txt                  legend with the numbers from this run
## ========================================================================= ##


# ---- 1. Settings ----------------------------------------------------------

base_dir     <- "C:/Users/mrzsjs1/Desktop/Addiction"
analysis_dir <- file.path(base_dir, "analysis")
out_dir      <- file.path(analysis_dir, "figures", "FigS1")
PHOS_DIR     <- analysis_dir                        # *_phosphoproteomics_significant.xlsx
METAB_XLSX   <- "C:/Users/mrzsjs1/Desktop/Addiction/analysis/metabolomics final.xlsx"
METAB_SHEET  <- "Filtered"
DET_DIR      <- analysis_dir                        # annotated DESeq2 result files

## kallisto output folders; each sample's .h5 file is found by searching these
## folders and reading the output name recorded inside the file ("-o .../W1"),
## so the location of each file within the folder does not matter
KNS_DIR <- "C:/Users/mrzsjs1/Desktop/Addiction/KNS42data"
SVG_DIR <- "C:/Users/mrzsjs1/Desktop/Addiction/SVGAdata"

## media-only control samples: sample ID as recorded by kallisto, and cell line
SAMPLE_IDS <- tibble::tribble(
  ~sample, ~line,    ~dir,
  "W1",    "KNS-42", KNS_DIR,
  "W2",    "KNS-42", KNS_DIR,
  "W3",    "KNS-42", KNS_DIR,
  "UC1",   "SVG-A",  SVG_DIR,
  "UC2",   "SVG-A",  SVG_DIR,
  "UC3",   "SVG-A",  SVG_DIR
)
## to set files by hand instead, give full paths here (same order as SAMPLE_IDS)
SAMPLE_FILES <- NULL

## optional pre-built transcript -> gene map (columns TXNAME, GENEID); used only
## if EnsDb.Hsapiens.v86 is not installed
TX2GENE_CSV <- NA


LINES <- c("KNS-42" = "KNS", "SVG-A" = "SVG")
LINE_ORDER <- c("KNS-42", "SVG-A")
DRUGS <- c(ethanol = "ethanol", THC = "thc", cocaine = "cocaine",
           morphine = "morphine", psilocin = "psilocin")
SVG_DRUGS <- c("ethanol", "THC", "cocaine", "psilocin")   # SVG-A morphine excluded
USE_CORRECTED <- TRUE
tx_file <- function(prefix, drug) {
  if (prefix == "KNS" && drug %in% c("ethanol", "psilocin"))
    return(sprintf("KNS_%s_DESeq2_transcript_results_CORRECTED.csv", drug))
  sprintf("%s_%s_DESeq2_transcript_results.csv", prefix, drug)
}
phos_file <- function(drug) sprintf("%s_phosphoproteomics_significant.xlsx", DRUGS[[drug]])
## drug group vs its control (THC vs ethanol vehicle); the sheet reports
## control / drug ratios, so log2FC is negated to give drug vs control
METAB_CONTRAST <- list(ethanol = c("E","UC"), THC = c("T","EC"), cocaine = c("C","UC"),
                       morphine = c("M","UC"), psilocin = c("P","UC"))

PADJ <- 0.05; LFC <- 1
NEAR_P <- 0.05; NEAR_LFC <- 0.5             # near-threshold definition (panel C context)
thr_lab <- sprintf("padj < %.2g & |log2FC| > %g", PADJ, LFC)
ROBUST_TPM <- 1; LOW_COUNTS <- 10
PADJ <- 0.05; LFC <- 1                     # drug-regulation rule (as in the paper)
USE_CORRECTED <- TRUE

LFC_CAP   <- 3        # colour scale limit, gene-level values
SHIFT_CAP <- 1        # colour scale limit, panel B

TARGETS <- tibble::tribble(
  ~gene,     ~class,
  "CNR1",    "Cannabinoid (THC)",
  "CNR2",    "Cannabinoid (THC)",
  "GPR55",   "Cannabinoid (THC)",
  "TRPV1",   "Cannabinoid (THC)",
  "OPRM1",   "Opioid (morphine)",
  "OPRD1",   "Opioid (morphine)",
  "OPRK1",   "Opioid (morphine)",
  "HTR1A",   "Serotonin (psilocin)",
  "HTR1B",   "Serotonin (psilocin)",
  "HTR1D",   "Serotonin (psilocin)",
  "HTR1E",   "Serotonin (psilocin)",
  "HTR1F",   "Serotonin (psilocin)",
  "HTR2A",   "Serotonin (psilocin)",
  "HTR2B",   "Serotonin (psilocin)",
  "HTR2C",   "Serotonin (psilocin)",
  "SLC6A3",  "Cocaine targets",
  "SLC6A2",  "Cocaine targets",
  "SLC6A4",  "Cocaine targets",
  "SIGMAR1", "Cocaine targets",
  "GRIN1",   "NMDA / GABA-A (ethanol)",
  "GRIN2A",  "NMDA / GABA-A (ethanol)",
  "GRIN2B",  "NMDA / GABA-A (ethanol)",
  "GABRA1",  "NMDA / GABA-A (ethanol)",
  "GABRB3",  "NMDA / GABA-A (ethanol)",
  "SLC1A2",  "Astrocytic glutamate transport",
  "SLC1A3",  "Astrocytic glutamate transport",
  "GFAP",    "Reference",
  "GAPDH",   "Reference"
)


## hallmark sets. Genes are the canonical readouts of each response; edit freely.
SETS <- list(
  "p53 / DNA damage" = c("CDKN1A","MDM2","GADD45A","BBC3","PMAIP1","TP53I3","SESN1","SESN2",
                         "RRM2B","TNFRSF10B","BAX","FAS","TP53INP1","ZMAT3","TRIAP1","DDB2",
                         "XPC","AEN","PLK3","TIGAR"),
  "Apoptosis effectors" = c("CASP3","CASP7","CASP8","CASP9","APAF1","BID","BAK1","BCL2L11",
                            "DIABLO","CYCS"),
  "ER stress / UPR" = c("DDIT3","HSPA5","ATF3","ATF4","XBP1","ERN1","EIF2AK3","PPP1R15A",
                        "TRIB3","HERPUD1","DNAJB9","EDEM1","ASNS","CHAC1"),
  "Oxidative stress (NRF2)" = c("HMOX1","NQO1","GCLM","GCLC","TXNRD1","SLC7A11","SRXN1","G6PD",
                                "OSGIN1","FTL","FTH1","PRDX1","GPX1","SOD1","SOD2","CAT",
                                "MT1X","MT2A"),
  "Heat shock" = c("HSPA1A","HSPA1B","HSPH1","DNAJA1","DNAJB1","HSPB1","HSP90AA1","BAG3"),
  "Proliferation / cell cycle" = c("MKI67","TOP2A","PCNA","CCNB1","CCNA2","CDK1","BUB1","AURKA",
                                   "AURKB","PLK1","RRM2","TYMS","MCM2","MCM3","MCM5","MCM7",
                                   "E2F1","FOXM1","KIF11","CENPF","TK1","CDC20","UBE2C","BIRC5")
)


# text sizes in printed points (figure drawn 180 mm wide)
FIG_W_MM <- 180; TEXT_SCALE <- 1
PT <- function(x) x * TEXT_SCALE
SZ <- list(tick = 7, axis = 8, strip = 7.5, title = 8.5, legend = 7, annot = 6.5, tag = 12)
fig_w   <- FIG_W_MM / 25.4
HEIGHTS <- c(A = 5.0, B = 1.75, C = 2.4, DE = 2.9)        # panel heights, inches
PANEL_CAPTIONS <- TRUE      # print each panel legend beneath its panel
WRITE_XLSX     <- TRUE      # supplementary tables as one labelled workbook
WRITE_CSV      <- FALSE     # also write each table as a separate csv
CAPTION_CHARS  <- 132       # wrap width for the panel captions
CAPTION_LINE   <- 0.165     # height allowed per caption line, inches
LINE_COLS <- c("KNS-42" = "#1D9E75", "SVG-A" = "#7F77DD")
COL_UP <- "#E24B4A"; COL_DOWN <- "#378ADD"; COL_NS <- "#B4B2A9"
COL_CONC <- "#1D9E75"; COL_DISC <- "#E24B4A"


# ---- 2. Packages ----------------------------------------------------------

cran <- c("readr","readxl","dplyr","tidyr","tibble","ggplot2","scales","patchwork","ggrepel","openxlsx")
missing <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing)
for (p in c("tximport","DESeq2"))
  if (!requireNamespace(p, quietly = TRUE)) stop("Install Bioconductor package: ", p)
suppressPackageStartupMessages({
  library(readr); library(readxl); library(dplyr); library(tidyr); library(tibble)
  library(ggplot2); library(scales); library(patchwork); library(ggrepel)
  library(openxlsx); library(tximport); library(DESeq2)
})
select <- dplyr::select; filter <- dplyr::filter; rename <- dplyr::rename
mutate <- dplyr::mutate; count <- dplyr::count; slice <- dplyr::slice
slice_min <- dplyr::slice_min; desc <- dplyr::desc
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
if (!file.exists(METAB_XLSX)) stop("Metabolomics file not found: ", METAB_XLSX)
strip_ver <- function(x) sub("\\.[0-9]+$", "", x)


# ---- 3. Theme --------------------------------------------------------------

theme_pub <- theme_classic(base_size = PT(SZ$axis)) +
  theme(plot.title   = element_text(face = "bold", hjust = 0.5, size = PT(SZ$title),
                                    margin = margin(b = 3)),
        axis.title   = element_text(size = PT(SZ$axis)),
        axis.text    = element_text(colour = "grey10", size = PT(SZ$tick)),
        strip.text   = element_text(face = "bold", size = PT(SZ$strip), margin = margin(2, 2, 2, 2)),
        strip.background = element_rect(fill = "grey95", colour = NA),
        plot.caption = element_text(size = PT(SZ$annot), colour = "grey40", hjust = 0),
        legend.position = "bottom", legend.title = element_text(size = PT(SZ$legend), vjust = 0.9),
        legend.text  = element_text(size = PT(SZ$legend)),
        legend.key.height = unit(6, "pt"), legend.key.width = unit(18, "pt"),
        legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(3, "pt"),
        plot.title.position = "panel", plot.caption.position = "plot",
        plot.margin  = margin(4, 6, 4, 4))
title_theme <- theme(plot.title = element_text(face = "bold", hjust = 0.5, size = PT(SZ$title)))

## a panel legend, rendered as its own row beneath the panel
cap_wrap <- function(letter, text)
  paste(strwrap(paste0("(", letter, ") ", text), CAPTION_CHARS), collapse = "\n")
cap_plot <- function(txt) {
  ggplot() +
    annotate("text", x = 0, y = 1, label = txt, hjust = 0, vjust = 1,
             size = PT(SZ$annot) / .pt, colour = "grey25", lineheight = 1.2) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme_void() + theme(plot.margin = margin(1, 4, 4, 4))
}
cap_rows <- function(txt) length(strsplit(txt, "\n")[[1]])

GEOM_TEXT <- PT(SZ$annot) / .pt
## supplementary tables are collected here and written once, at the end
SUPP <- list()
add_table <- function(sheet, number, title, description, data) {
  SUPP[[sheet]] <<- list(number = number, title = title,
                         description = description, data = as.data.frame(data))
  if (WRITE_CSV) write_csv(data, file.path(out_dir, paste0(sheet, ".csv")))
  invisible(data)
}

write_supplementary <- function(file) {
  wb <- createWorkbook()
  hdr  <- createStyle(textDecoration = "bold", fgFill = "#E8E8E8", border = "bottom",
                      halign = "left", valign = "center", wrapText = FALSE)
  ttl  <- createStyle(textDecoration = "bold", fontSize = 12)
  note <- createStyle(fontColour = "#404040", wrapText = FALSE)

  ## contents sheet
  addWorksheet(wb, "Contents")
  writeData(wb, "Contents", "Supplementary Tables", startCol = 1, startRow = 1)
  addStyle(wb, "Contents", ttl, rows = 1, cols = 1)
  idx <- data.frame(Sheet = names(SUPP),
                    Table = vapply(SUPP, `[[`, "", "number"),
                    Title = vapply(SUPP, `[[`, "", "title"),
                    Description = vapply(SUPP, `[[`, "", "description"),
                    Rows = vapply(SUPP, function(x) nrow(x$data), 1L))
  writeData(wb, "Contents", idx, startRow = 3, headerStyle = hdr)
  setColWidths(wb, "Contents", cols = 1:5, widths = c(26, 10, 58, 110, 8))
  freezePane(wb, "Contents", firstActiveRow = 4)

  ## one sheet per table: title, description, then the data
  for (nm in names(SUPP)) {
    x <- SUPP[[nm]]
    addWorksheet(wb, nm)
    writeData(wb, nm, paste0(x$number, ". ", x$title), startRow = 1)
    addStyle(wb, nm, ttl, rows = 1, cols = 1)
    writeData(wb, nm, x$description, startRow = 2)
    addStyle(wb, nm, note, rows = 2, cols = 1)
    writeData(wb, nm, x$data, startRow = 4, headerStyle = hdr)
    freezePane(wb, nm, firstActiveRow = 5)
    setColWidths(wb, nm, cols = seq_len(ncol(x$data)), widths = "auto")
    addFilter(wb, nm, rows = 4, cols = seq_len(ncol(x$data)))
  }
  saveWorkbook(wb, file, overwrite = TRUE)
  message(sprintf("Supplementary tables written: %s (%d tables)", file, length(SUPP)))
}

save_fig <- function(p, name, w, h) {
  ggsave(file.path(out_dir, paste0(name, ".pdf")), p, width = w, height = h,
         device = cairo_pdf, limitsize = FALSE)
  ggsave(file.path(out_dir, paste0(name, ".png")), p, width = w, height = h,
         dpi = 600, limitsize = FALSE, bg = "white")
}


# ---- 4. Transcript -> gene map --------------------------------------------

if (requireNamespace("EnsDb.Hsapiens.v86", quietly = TRUE)) {
  tx <- ensembldb::transcripts(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86,
                               columns = c("tx_id","gene_name"), return.type = "data.frame")
  tx2gene <- data.frame(TXNAME = tx$tx_id, GENEID = toupper(tx$gene_name))
  map_src <- "EnsDb.Hsapiens.v86"
} else if (!is.na(TX2GENE_CSV) && file.exists(TX2GENE_CSV)) {
  tx2gene <- read.csv(TX2GENE_CSV, stringsAsFactors = FALSE)[, c("TXNAME","GENEID")]
  map_src <- basename(TX2GENE_CSV)
} else {
  ## fallback: tested transcripts from the annotated DESeq2 files + org.Hs.eg.db
  fs <- list.files(DET_DIR, pattern = "DESeq2_transcript_results.*\\.csv$", full.names = TRUE)
  fs <- fs[!grepl("noGeneSymbol", fs)]
  m1 <- bind_rows(lapply(fs, function(f)
    read_csv(f, col_select = c("ENSEMBL_TRANSCRIPT","gene_symbol"), show_col_types = FALSE))) %>%
    transmute(TXNAME = strip_ver(ENSEMBL_TRANSCRIPT), GENEID = toupper(gene_symbol))
  m2 <- NULL
  if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    k  <- AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype = "ENSEMBLTRANS")
    m2 <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = k,
                                                 keytype = "ENSEMBLTRANS", columns = "SYMBOL")) %>%
      transmute(TXNAME = ENSEMBLTRANS, GENEID = toupper(SYMBOL))
  }
  tx2gene <- bind_rows(m1, m2)
  map_src <- "DESeq2 files + org.Hs.eg.db (partial; install EnsDb.Hsapiens.v86 for full coverage)"
}
tx2gene <- tx2gene %>% filter(!is.na(GENEID), GENEID != "") %>%
  mutate(TXNAME = strip_ver(TXNAME)) %>% distinct(TXNAME, .keep_all = TRUE) %>% as.data.frame()
message("Transcript -> gene map: ", map_src, " (", nrow(tx2gene), " transcripts)")
unmapped_targets <- setdiff(TARGETS$gene, tx2gene$GENEID)
if (length(unmapped_targets))
  message("Targets absent from the map (reported as 'n.a.'): ", paste(unmapped_targets, collapse = ", "))


# ---- 5. Import ------------------------------------------------------------

if (!requireNamespace("rhdf5", quietly = TRUE))
  stop("Reading kallisto .h5 files needs rhdf5: BiocManager::install('rhdf5')")

## sample ID recorded in a kallisto .h5 file (the name of its -o output folder)
h5_sample <- function(f) {
  call <- tryCatch(as.character(rhdf5::h5read(f, "aux/call"))[1], error = function(e) NA_character_)
  rhdf5::h5closeAll()
  if (is.na(call)) return(NA_character_)
  toks <- strsplit(call, "\\s+")[[1]]
  i <- which(toks == "-o")[1]
  if (is.na(i)) return(NA_character_)
  basename(toks[i + 1])
}
locate_sample <- function(id, dir) {
  if (!dir.exists(dir)) stop("Folder not found: ", dir)
  h5 <- list.files(dir, pattern = "\\.h5$", recursive = TRUE, full.names = TRUE)
  ## check likely files first (folder or file named after the sample)
  h5 <- h5[order(!grepl(paste0("(^|/)", id, "([_/.]|$)"), h5, ignore.case = TRUE))]
  for (f in h5) {
    o <- h5_sample(f)
    if (!is.na(o) && grepl(paste0("^", id, "($|_)"), o)) return(f)
  }
  stop("No kallisto .h5 file for sample ", id, " found under ", dir)
}

SAMPLES <- SAMPLE_IDS %>%
  mutate(file = if (is.null(SAMPLE_FILES))
                  mapply(locate_sample, sample, dir, USE.NAMES = FALSE)
                else SAMPLE_FILES) %>%
  select(sample, line, file)
missing_files <- SAMPLES$file[!file.exists(SAMPLES$file)]
if (length(missing_files)) stop("Sample files not found:\n  ", paste(missing_files, collapse = "\n  "))
message("Samples used:")
for (i in seq_len(nrow(SAMPLES)))
  message(sprintf("  %-4s %-7s %s", SAMPLES$sample[i], SAMPLES$line[i], SAMPLES$file[i]))
write_csv(SAMPLES, file.path(out_dir, "baseline_samples_used.csv"))

SAMPLES <- SAMPLES %>%
  mutate(line = factor(line, levels = LINE_ORDER))
print(SAMPLES)
if (any(table(SAMPLES$line) < 2)) stop("Need >= 2 samples per line - check the file paths.")

## tximport only reads a kallisto file as HDF5 when it is named abundance.h5,
## so each file is copied to a temporary folder under that name
stage_dir <- file.path(tempdir(), "baseline_kallisto")
files <- setNames(vapply(seq_len(nrow(SAMPLES)), function(i) {
  f <- SAMPLES$file[i]
  if (!grepl("\\.h5$", f) || basename(f) == "abundance.h5") return(f)
  d <- file.path(stage_dir, SAMPLES$sample[i])
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  file.copy(f, file.path(d, "abundance.h5"), overwrite = TRUE)
  file.path(d, "abundance.h5")
}, character(1)), SAMPLES$sample)
txi <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreTxVersion = TRUE)
TPM <- txi$abundance; CNT <- txi$counts
write_csv(as.data.frame(TPM) %>% rownames_to_column("gene"), file.path(out_dir, "baseline_gene_TPM.csv"))
write_csv(as.data.frame(round(CNT)) %>% rownames_to_column("gene"), file.path(out_dir, "baseline_gene_counts.csv"))
mapped_pct <- round(100 * colSums(TPM), 1) / 1e6
message("TPM assigned to a gene: ", paste(sprintf("%s %.1f%%", names(mapped_pct), mapped_pct), collapse = ", "))


# ---- 6. Per-line summaries -------------------------------------------------

get <- function(M, g) if (g %in% rownames(M)) M[g, ] else setNames(rep(0, ncol(M)), colnames(M))
long <- bind_rows(lapply(TARGETS$gene, function(g)
  tibble(gene = g, sample = colnames(TPM), TPM = get(TPM, g), counts = get(CNT, g)))) %>%
  left_join(SAMPLES %>% select(sample, line), by = "sample")

pct <- bind_rows(lapply(LINE_ORDER, function(ln) {
  m <- rowMeans(TPM[, SAMPLES$sample[SAMPLES$line == ln], drop = FALSE])
  m <- m[m > 0]
  tibble(gene = names(m), line = ln, percentile = round(100 * rank(m) / length(m), 1))
}))

summ <- long %>% group_by(gene, line) %>%
  summarise(mean_TPM = mean(TPM), sd_TPM = sd(TPM), n = n(),
            min_TPM = min(TPM), min_counts = min(counts),
            counts = paste(round(counts), collapse = "/"), .groups = "drop") %>%
  mutate(line = as.character(line),
         tier = case_when(gene %in% unmapped_targets ~ "n.a.",
                          min_TPM >= ROBUST_TPM ~ "robust",
                          min_counts >= LOW_COUNTS ~ "low",
                          TRUE ~ "absent")) %>%
  left_join(pct, by = c("gene","line"))


# ---- 7. KNS-42 vs SVG-A (DESeq2 on the control samples) -------------------

cd <- data.frame(line = relevel(droplevels(SAMPLES$line), ref = "SVG-A"), row.names = SAMPLES$sample)
dds <- DESeqDataSetFromTximport(txi, colData = cd, design = ~ line)
dds <- dds[rowSums(counts(dds)) > 0, ]
dds <- DESeq(dds, quiet = TRUE)
rn  <- grep("^line_", resultsNames(dds), value = TRUE)[1]
res <- results(dds, name = rn)                           # BH across all genes
line_stats <- as.data.frame(res) %>% rownames_to_column("gene") %>%
  transmute(gene, log2FC_KNS_vs_SVG = round(log2FoldChange, 2),
            lfcSE = round(lfcSE, 2), pvalue = signif(pvalue, 3), padj = signif(padj, 3))


# ---- 8. Drug regulation of each target (existing DESeq2 results) ---------

find_det <- function(prefix, drug) {
  f <- list.files(DET_DIR, pattern = "\\.csv$", full.names = TRUE)
  f <- f[!grepl("noGeneSymbol", basename(f), ignore.case = TRUE)]
  hit <- f[grepl(paste0("^", prefix, "_"), basename(f), ignore.case = TRUE) &
             grepl(paste0("_", drug, "_"), basename(f), ignore.case = TRUE)]
  if (!length(hit)) return(NA_character_)
  corr <- grepl("CORRECTED", basename(hit), ignore.case = TRUE)
  if (USE_CORRECTED && any(corr)) hit[corr][1] else hit[!corr][1]
}
reg_drugs <- names(DRUGS)
exps <- bind_rows(lapply(c("KNS-42" = "KNS", "SVG-A" = "SVG"), function(pre)
  tibble(experiment = reg_drugs,
         file = vapply(reg_drugs, function(d) find_det(pre, DRUGS[[d]]), ""))), .id = "line") %>%
  filter(!is.na(file), !(line == "SVG-A" & experiment == "morphine"))   # excluded in the paper
regw <- bind_rows(lapply(seq_len(nrow(exps)), function(i)
  read_csv(exps$file[i], show_col_types = FALSE,
           col_select = c("log2FoldChange","padj","gene_symbol")) %>%
    mutate(gene = toupper(gene_symbol)) %>%
    filter(gene %in% TARGETS$gene, !is.na(padj)) %>%
    group_by(gene) %>% slice_min(padj, n = 1, with_ties = FALSE) %>% ungroup() %>%
    transmute(gene, line = exps$line[i], experiment = exps$experiment[i],
              log2FC = log2FoldChange, padj, sig = padj < PADJ & abs(log2FC) > LFC)))
reg_txt <- regw %>% filter(sig) %>%
  mutate(txt = sprintf("%s %s (%+.1f, padj %.2g)", line, experiment, log2FC, padj)) %>%
  group_by(gene) %>% summarise(regulated_by = paste(txt, collapse = "; "), .groups = "drop")


# ---- 9. Table and text summary --------------------------------------------

k <- function(ln) ifelse(ln == "KNS-42", "KNS42", "SVGA")
wide <- summ %>% mutate(line = k(line)) %>%
  pivot_wider(id_cols = gene, names_from = line,
              values_from = c(mean_TPM, sd_TPM, counts, tier, percentile),
              names_glue = "{line}_{.value}")
tab <- TARGETS %>% rename(target_class = class) %>%
  left_join(wide, by = "gene") %>% left_join(line_stats, by = "gene") %>%
  left_join(reg_txt, by = "gene") %>%
  mutate(across(matches("TPM$"), ~ round(.x, 2)),
         regulated_by = coalesce(regulated_by, ""),
         across(c(log2FC_KNS_vs_SVG, padj), ~ ifelse(gene %in% unmapped_targets, NA, .x)))
add_table("S1_baseline_targets", "Table S1",
  "Baseline expression of canonical drug targets",
  paste("Media-only control cells (n = 3 per line). Mean +/- SD TPM, read counts per replicate, expression call",
        "(robust: TPM >= 1 in all replicates; low: TPM < 1 but >= 10 reads in all replicates; absent: neither)",
        "and within-line percentile; DESeq2 log2 fold-change and adjusted p for KNS-42 versus SVG-A; and the drug",
        "treatments that regulated each target (adjusted p < 0.05, |log2FC| > 1)."),
  tab)

f1 <- function(m, s) ifelse(is.na(m), "NA", sprintf("%.2f +/- %.2f", m, s))
txt <- with(tab, sprintf("%-8s KNS-42 %s TPM [%s; reads %s] | SVG-A %s TPM [%s; reads %s] | log2FC %s, padj %s%s",
  gene, f1(KNS42_mean_TPM, KNS42_sd_TPM), KNS42_tier, KNS42_counts,
  f1(SVGA_mean_TPM, SVGA_sd_TPM), SVGA_tier, SVGA_counts,
  ifelse(is.na(log2FC_KNS_vs_SVG), "NA", sprintf("%+.1f", log2FC_KNS_vs_SVG)),
  ifelse(is.na(padj), "NA", format(padj, digits = 2)),
  ifelse(regulated_by == "", "", paste0(" | regulated: ", regulated_by))))
writeLines(c(sprintf("Media-only controls: %s. Transcript-to-gene map: %s.",
                     paste(sprintf("%s n = %d", names(table(SAMPLES$line)), table(SAMPLES$line)), collapse = ", "),
                     map_src),
             sprintf("Tiers: robust = TPM >= %g in all samples; low = TPM < %g but >= %d reads in all samples.",
                     ROBUST_TPM, ROBUST_TPM, LOW_COUNTS),
             txt), file.path(out_dir, "baseline_text_summary.txt"))
cat(txt, sep = "\n")


# ---- 10. Panel A -----------------------------------------------------------

lev_gene  <- rev(TARGETS$gene)
lev_class <- unique(TARGETS$class)
theme_base <- theme_classic(base_size = PT(SZ$axis)) +
  theme(axis.text   = element_text(colour = "grey10", size = PT(SZ$tick)),
        axis.text.y = element_text(face = "italic"),
        strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold", size = PT(SZ$strip)),
        strip.background = element_blank(), strip.placement = "outside",
        panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.25),
        panel.spacing.y = unit(3, "pt"),
        legend.position = "top", legend.title = element_blank(),
        legend.text = element_text(size = PT(SZ$legend)),
        legend.margin = margin(0, 0, 0, 0), plot.margin = margin(4, 2, 4, 4))
facet_rows <- function(...) facet_grid(class ~ ., scales = "free_y", space = "free_y", ...)

pd <- long %>% filter(!gene %in% unmapped_targets) %>% left_join(TARGETS, by = "gene") %>%
  mutate(gene = factor(gene, levels = lev_gene), class = factor(class, levels = lev_class),
         y = log2(TPM + 1))
md <- pd %>% group_by(gene, class, line) %>% summarise(y = mean(y), .groups = "drop")
blank_rows <- TARGETS %>% mutate(gene = factor(gene, levels = lev_gene),
                                 class = factor(class, levels = lev_class))
dodge <- position_dodge(width = 0.6)

p_main <- ggplot(pd, aes(y, gene, colour = line)) +
  geom_blank(data = blank_rows, aes(x = 0, y = gene), inherit.aes = FALSE) +
  geom_vline(xintercept = log2(ROBUST_TPM + 1), linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  geom_crossbar(data = md, aes(xmin = y, xmax = y), width = 0.55, linewidth = 0.35,
                position = dodge, show.legend = FALSE) +
  geom_point(position = position_jitterdodge(jitter.height = 0.1, jitter.width = 0,
                                             dodge.width = 0.6, seed = 1),
             size = 1.3, alpha = 0.85, stroke = 0) +
  scale_colour_manual(values = LINE_COLS) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.03)),
                     sec.axis = dup_axis(name = "TPM",
                                         breaks = log2(c(1, 10, 100, 1000) + 1),
                                         labels = c("1", "10", "100", "1000"))) +
  facet_rows(switch = "y") +
  guides(colour = guide_legend(override.aes = list(size = 2.5))) +
  labs(x = expression(log[2]*"(TPM + 1)"), y = NULL) +
  theme_base

sym <- function(t) recode(t, robust = "+", low = "(+)", absent = "\u2013", n.a. = "n.a.")
sd_df <- tab %>% left_join(TARGETS, by = "gene") %>%
  transmute(gene = factor(gene, levels = lev_gene), class = factor(class, levels = lev_class),
            det = paste0(sym(KNS42_tier), " / ", sym(SVGA_tier)),
            none = KNS42_tier %in% c("absent","n.a.") & SVGA_tier %in% c("absent","n.a."),
            ## no fold change is shown when the gene is not detected in either line
            fc  = ifelse(is.na(log2FC_KNS_vs_SVG) | none, "", sprintf("%+.1f", log2FC_KNS_vs_SVG)),
            pl  = case_when(gene %in% unmapped_targets | none ~ "",
                            is.na(padj) ~ "n.t.", padj < 0.001 ~ "***",
                            padj < 0.01 ~ "**", padj < 0.05 ~ "*", TRUE ~ "ns"))
p_stat <- ggplot(sd_df, aes(y = gene)) +
  geom_text(aes(x = 0.75, label = det), size = PT(SZ$annot) / .pt) +
  geom_text(aes(x = 2.6,  label = fc),  size = PT(SZ$annot) / .pt) +
  geom_text(aes(x = 4.0,  label = pl),  size = PT(SZ$annot) / .pt) +
  scale_x_continuous(limits = c(-0.3, 4.6), breaks = c(0.75, 2.6, 4.0), position = "top",
                     labels = c("Expressed\nK / S", "log2FC\nK vs S", "padj")) +
  facet_rows() +
  theme_void() +
  theme(axis.text.x.top = element_text(size = PT(SZ$annot), lineheight = 0.9, vjust = 0),
        strip.text = element_blank(), panel.spacing.y = unit(3, "pt"),
        plot.margin = margin(4, 2, 4, 0))

rg <- expand_grid(TARGETS, exps %>% select(line, experiment)) %>%
  left_join(regw, by = c("gene","line","experiment")) %>%
  mutate(col = paste0(ifelse(line == "KNS-42", "K ", "S "),
                      recode(experiment, ethanol = "EtOH", cocaine = "Coc", morphine = "Mor", psilocin = "Psi")),
         fill = ifelse(!is.na(sig) & sig, pmax(pmin(log2FC, 4), -4), NA),
         gene = factor(gene, levels = lev_gene), class = factor(class, levels = lev_class))
rg$col <- factor(rg$col, levels = unique(rg$col[order(match(rg$line, LINE_ORDER), match(rg$experiment, names(DRUGS)))]))
p_reg <- ggplot(rg, aes(col, gene, fill = fill)) +
  geom_tile(colour = "grey85", linewidth = 0.25) +
  scale_fill_gradient2(low = "#378ADD", mid = "white", high = "#E24B4A", midpoint = 0,
                       limits = c(-4, 4), na.value = "white", name = "Drug\nlog2FC",
                       breaks = c(-4, 0, 4)) +
  scale_x_discrete(position = "top") +
  facet_rows() + labs(x = NULL, y = NULL) +
  theme_minimal(base_size = PT(SZ$axis)) +
  theme(axis.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5, size = PT(SZ$annot)),
        axis.text.y = element_blank(), panel.grid = element_blank(),
        strip.text = element_blank(), panel.spacing.y = unit(3, "pt"),
        legend.position = "right", legend.title = element_text(size = PT(SZ$annot)),
        legend.text = element_text(size = PT(SZ$annot)),
        legend.key.height = unit(14, "pt"), legend.key.width = unit(6, "pt"),
        plot.margin = margin(4, 4, 4, 2))

cap_A <- paste0(
  "Expression of the canonical molecular targets of each substance in media-only control cells (n = 3 per line), shown as log2(TPM + 1). ",
  "Points are individual replicates and bars are means. The dashed line marks TPM = 1. A gene is marked + where TPM was at least 1 in every ",
  "replicate, (+) where TPM was below 1 but at least 10 reads were present in every replicate, and - where it was not detected ",
  "(K, KNS-42; S, SVG-A). Log2 fold-change and adjusted p compare KNS-42 with SVG-A (DESeq2, Benjamini-Hochberg; * p < 0.05, ** p < 0.01, ",
  "*** p < 0.001). The two cell lines were sequenced in separate runs, so these comparisons are descriptive. The grid on the right shows ",
  "which drugs regulated each target (adjusted p < 0.05, |log2FC| > 1).")
panel_A <- p_main + p_stat + p_reg + plot_layout(widths = c(1.9, 1.75, 1.05))


# ---- 11. Read gene-level results --------------------------------------------
## (find_det is defined in section 8)

cyto_arms <- bind_rows(lapply(names(LINES), function(ln) {
  dr <- if (ln == "SVG-A") SVG_DRUGS else names(DRUGS)
  tibble(line = ln, drug = dr,
         file = vapply(dr, function(d) find_det(LINES[[ln]], DRUGS[[d]]), character(1)))
})) %>% filter(!is.na(file))
if (!nrow(cyto_arms)) stop("No DESeq2 result files found in ", analysis_dir)
print(as.data.frame(cyto_arms %>% mutate(file = basename(file))), row.names = FALSE)

## one value per gene: the most significant transcript
GENES <- bind_rows(lapply(seq_len(nrow(cyto_arms)), function(i) {
  read_csv(cyto_arms$file[i], show_col_types = FALSE,
           col_select = c("log2FoldChange","padj","gene_symbol")) %>%
    mutate(gene = toupper(trimws(gene_symbol))) %>%
    filter(!is.na(gene), gene != "", !is.na(padj), !is.na(log2FoldChange)) %>%
    arrange(padj) %>% distinct(gene, .keep_all = TRUE) %>%
    transmute(line = cyto_arms$line[i], drug = cyto_arms$drug[i], gene,
              lfc = log2FoldChange, padj)
}))

set_tbl <- tibble(set = rep(names(SETS), lengths(SETS)),
                  gene = unlist(SETS, use.names = FALSE)) %>%
  mutate(set = factor(set, levels = names(SETS)))


# ---- 12. Set-level tests ------------------------------------------------------

## competitive test: hallmark genes against every other quantified gene
tests <- bind_rows(lapply(seq_len(nrow(cyto_arms)), function(i) {
  g <- GENES %>% filter(line == cyto_arms$line[i], drug == cyto_arms$drug[i])
  bind_rows(lapply(names(SETS), function(s) {
    inset <- g$gene %in% SETS[[s]]
    a <- g$lfc[inset]; b <- g$lfc[!inset]
    sig <- g %>% filter(inset, padj < PADJ, abs(lfc) > LFC)
    tibble(line = cyto_arms$line[i], drug = cyto_arms$drug[i], set = s,
           n_tested = length(a), n_sig = nrow(sig),
           n_up = sum(sig$lfc > 0), n_down = sum(sig$lfc < 0),
           median_set = median(a), median_background = median(b),
           p = if (length(a) > 2) suppressWarnings(wilcox.test(a, b)$p.value) else NA_real_)
  }))
})) %>%
  mutate(padj_set = p.adjust(p, "BH"),
         set = factor(set, levels = names(SETS)),
         drug = factor(drug, levels = names(DRUGS)),
         line = factor(line, levels = names(LINES)))
add_table("S2a_cytotoxicity_sets", "Table S2a",
  "Cytotoxicity hallmark sets: set-level tests",
  paste("For every cell line and treatment: genes tested (n_tested), genes reaching adjusted p < 0.05 and",
        "|log2FC| > 1 (n_sig, n_up, n_down), the median log2 fold-change of the set (median_set) and of all other",
        "quantified genes (median_background), and the two-sided Wilcoxon rank-sum p-value (p) with its",
        "Benjamini-Hochberg adjustment across all set-by-treatment tests (padj_set)."),
  tests)

gene_tbl <- GENES %>% inner_join(set_tbl, by = "gene") %>%
  mutate(drug = factor(drug, levels = names(DRUGS)),
         line = factor(line, levels = names(LINES)),
         sig = padj < PADJ & abs(lfc) > LFC) %>%
  arrange(set, gene)
add_table("S2b_cytotoxicity_gene_sets", "Table S2b",
  "Cytotoxicity hallmark sets: membership",
  "Genes assigned to each hallmark set, and whether each was quantified in the transcriptomic data.",
  set_tbl %>% mutate(quantified = gene %in% GENES$gene))

add_table("S2c_cytotoxicity_genes", "Table S2c",
  "Cytotoxicity hallmark genes: per-treatment values",
  paste("Log2 fold-change (lfc) and adjusted p-value of every hallmark gene in every cell line and treatment;",
        "sig marks genes reaching adjusted p < 0.05 and |log2FC| > 1."),
  gene_tbl)
## the hallmark sets themselves, including genes not quantified in a given arm

n_tests <- sum(!is.na(tests$p)); n_hit <- sum(tests$padj_set < 0.05, na.rm = TRUE)
message(sprintf("\n%d of %d set-level tests significant at FDR < 0.05", n_hit, n_tests))
if (n_hit) print(as.data.frame(tests %>% filter(padj_set < 0.05) %>%
                                 mutate(across(where(is.numeric), ~ signif(.x, 3)))), row.names = FALSE)
message(sprintf("Hallmark genes passing %s: %d of %d gene x treatment tests",
                sprintf("padj < %.2g and |log2FC| > %g", PADJ, LFC),
                sum(gene_tbl$sig), nrow(gene_tbl)))


# ---- 13. Panel B -----------------------------------------------------------

pd_S <- tests %>%
  mutate(shift = median_set - median_background,
         arm = factor(as.character(drug), levels = names(DRUGS)),
         set = factor(set, levels = rev(names(SETS))),
         hit = !is.na(padj_set) & padj_set < 0.05)
p_SUM <- ggplot(pd_S, aes(arm, set, fill = pmax(pmin(shift, SHIFT_CAP), -SHIFT_CAP))) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_tile(data = filter(pd_S, hit), fill = NA, colour = "grey10", linewidth = 0.6) +
  geom_text(aes(label = n_sig), size = PT(SZ$annot) / .pt, colour = "grey15") +
  scale_fill_gradient2(low = COL_DOWN, mid = "white", high = COL_UP, midpoint = 0,
                       limits = c(-SHIFT_CAP, SHIFT_CAP),
                       name = expression("Set shift ("*log[2]*"FC)")) +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
  facet_grid(~ line, scales = "free_x", space = "free_x") +
  labs(title = "Cytotoxicity and stress hallmark genes at 24 h", x = NULL, y = NULL,
       caption = sprintf("Colour: median log2FC of the set minus all other genes; numbers: set genes with adjusted p < %.2g and |log2FC| > %g;\noutline: FDR < 0.05 (Wilcoxon rank-sum, Benjamini-Hochberg; %d of %d tests).",
                         PADJ, LFC, n_hit, n_tests)) +
  theme_pub +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.y = element_text(size = PT(SZ$annot)),
        panel.spacing.x = unit(6, "pt"))
cap_B <- paste0(
  "Hallmark transcriptional responses to cytotoxic exposure, tested in each treatment: p53/DNA-damage signalling, apoptosis effectors, the ",
  "unfolded protein response, NRF2-dependent oxidative stress, heat-shock genes, and proliferation and cell-cycle genes. Colour shows the ",
  "median log2 fold-change of each set minus that of all other quantified genes. Numbers give the set genes reaching adjusted p < 0.05 and ",
  "|log2FC| > 1. Outlined cells differ from background at FDR < 0.05 (two-sided Wilcoxon rank-sum, Benjamini-Hochberg; ",
  sprintf("%d of %d tests). Gene-level values are given in Supplementary Table S2.", n_hit, n_tests))
panel_B <- p_SUM + labs(caption = NULL)   # the panel legend is printed below


# ---- 14. Drug-unique features ----------------------------------------------

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


## drug-unique transcripts and genes

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


## phosphoproteomics and metabolomics (KNS-42)

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


# ---- 15. Supplementary tables for panels C-E --------------------------------
## gene, phosphoprotein and metabolite lists behind panels C-E

add_table("S3a_drug_unique_counts", "Table S3a",
  "Drug-unique and shared features",
  paste("Numbers of features regulated by each drug within a cell line and omic layer (total), those regulated by",
        "that drug only (unique) and those also regulated by another drug (shared). SVG-A morphine was excluded."),
  counts)
add_table("S3b_drug_unique_genes", "Table S3b",
  "Drug-unique genes",
  paste("Genes significant for one drug only (adjusted p < 0.05, |log2FC| > 1) within a cell line.",
        "near_threshold marks genes that also showed a sub-threshold change (unadjusted p < 0.05, |log2FC| > 0.5)",
        "with at least one other drug."),
  bind_rows(lapply(names(UNIQ), function(ln)
    bind_rows(lapply(names(UNIQ[[ln]]), function(d)
      tibble(line = ln, drug = d, gene = UNIQ[[ln]][[d]]))))) %>%
      left_join(near_tab, by = c("line", "drug", "gene")))

add_table("S3c_drug_unique_phosphoproteins", "Table S3c",
  "Drug-unique phosphoproteins",
  "Phosphoproteins with at least one regulated phosphopeptide (unadjusted p < 0.05) for one drug only, in KNS-42.",
  bind_rows(lapply(names(uP), function(d)
    tibble(line = "KNS-42", drug = d, phosphoprotein = uP[[d]]))))

add_table("S3d_drug_unique_metabolites", "Table S3d",
  "Drug-unique metabolites",
  paste("Metabolite features significantly changed by one drug only (adjusted p < 0.05, |log2FC| > 1) in KNS-42.",
        "Positive log2FC indicates an increase with treatment."),
  bind_rows(lapply(names(DRUGS), function(d)
    METAB[[d]] %>% filter(sig, metabolite %in% uM[[d]]) %>%
      transmute(line = "KNS-42", drug = d, metabolite, log2FC, padj))))

message(sprintf("Drug-unique features written: %d genes, %d phosphoproteins, %d metabolites",
                sum(lengths(unlist(UNIQ, recursive = FALSE))), length(unlist(uP)), length(unlist(uM))))


# ---- 15. Panel C -----------------------------------------------------------
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
p_HEAT <- make_heat(NULL)
cap_C <- paste0(
  "Log2 fold-change of the drug-unique genes (columns) in response to every drug (rows), grouped by the drug to which they are unique, for ",
  "KNS-42 (upper) and SVG-A (lower). Values are capped at plus or minus 5; grey indicates a gene that was not quantified for that drug. ",
  "Gene lists are given in Supplementary Table S3.")
panel_C <- p_HEAT


# ---- 16. Panel D -----------------------------------------------------------

## genes unique to psilocin in BOTH cell lines
psi_genes <- intersect(UNIQ[["KNS-42"]][["psilocin"]], UNIQ[["SVG-A"]][["psilocin"]])
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
                      labels = dir_lab, name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 2.5))) +
  labs(title = sprintf("Psilocin-unique genes in both lines (n = %d)", nrow(pd_E)),
       x = expression(log[2]~"fold-change"~"(KNS-42)"),
       y = expression(log[2]~"fold-change"~"(SVG-A)"),
       caption = sprintf("%s. Gene log2FC from the most significant transcript.", thr_lab)) +
  theme_pub
cap_D <- sprintf(paste0(
  "The %d genes that were unique to psilocin in both cell lines, plotted as log2 fold-change in KNS-42 against SVG-A. Points are coloured by ",
  "whether the two cell lines changed in the same direction (concordant, green) or in opposite directions (discordant, red)."), nrow(pd_E))
panel_D <- p_E + labs(caption = NULL)


# ---- 17. Panel E -----------------------------------------------------------

pd_F <- METAB$psilocin %>%
  filter(sig, metabolite %in% uM$psilocin) %>%
  group_by(metabolite) %>% slice_min(padj, n = 1, with_ties = FALSE) %>% ungroup() %>%
  mutate(label = ifelse(nchar(metabolite) > 35, paste0(substr(metabolite, 1, 33), "..."), metabolite),
         label = reorder(label, log2FC),
         cls = ifelse(log2FC > 0, "Up", "Down"))

p_F <- ggplot(pd_F, aes(log2FC, label, fill = cls)) +
  geom_col(colour = "grey20", linewidth = 0.2, width = 0.74) +
  geom_vline(xintercept = 0, colour = "grey20", linewidth = 0.3) +
  scale_fill_manual(values = c(Up = COL_UP, Down = COL_DOWN), breaks = c("Up", "Down"), name = NULL) +
  labs(title = "Psilocin-unique metabolites",
       x = expression(log[2]~"fold-change"), y = NULL, caption = thr_lab) +
  theme_pub + theme(axis.text.y = element_text(size = PT(SZ$annot)),
                    plot.margin = margin(4, 26, 4, 4))
cap_E <- paste0(
  "Metabolite features that were significantly changed by psilocin and by no other drug (adjusted p < 0.05, |log2FC| > 1), ranked by log2 ",
  "fold-change. Red indicates an increase and blue a decrease with treatment.")
panel_E <- p_F + labs(caption = NULL)



# ---- 18. Assemble Supplementary Figure S1 ----------------------------------

caps <- list(A = cap_wrap("A", cap_A), B = cap_wrap("B", cap_B), C = cap_wrap("C", cap_C),
             DE = paste(cap_wrap("D", cap_D), cap_wrap("E", cap_E), sep = "\n"))
row_DE <- wrap_elements(full = panel_D) | wrap_elements(full = panel_E)

if (PANEL_CAPTIONS) {
  cap_h <- vapply(caps, function(x) cap_rows(x) * CAPTION_LINE, numeric(1))
  figS1 <- wrap_elements(full = panel_A) / cap_plot(caps$A) /
    wrap_elements(full = panel_B) / cap_plot(caps$B) /
    wrap_elements(full = panel_C) / cap_plot(caps$C) /
    row_DE / cap_plot(caps$DE) +
    plot_layout(heights = c(HEIGHTS[["A"]],  cap_h[["A"]],
                            HEIGHTS[["B"]],  cap_h[["B"]],
                            HEIGHTS[["C"]],  cap_h[["C"]],
                            HEIGHTS[["DE"]], cap_h[["DE"]])) +
    plot_annotation(tag_levels = list(c("A", "", "B", "", "C", "", "D", "E", ""))) &
    theme(plot.tag = element_text(face = "bold", size = PT(SZ$tag)))
  fig_h <- sum(HEIGHTS) + sum(cap_h)
} else {
  figS1 <- wrap_elements(full = panel_A) / wrap_elements(full = panel_B) /
    wrap_elements(full = panel_C) / row_DE +
    plot_layout(heights = HEIGHTS) +
    plot_annotation(tag_levels = list(c("A", "B", "C", "D", "E"))) &
    theme(plot.tag = element_text(face = "bold", size = PT(SZ$tag)))
  fig_h <- sum(HEIGHTS)
}
save_fig(figS1, "FigS1_supplementary", fig_w, fig_h)


# ---- 19. Legend -------------------------------------------------------------

nK <- sum(SAMPLES$line == "KNS-42"); nS <- sum(SAMPLES$line == "SVG-A")
legend <- paste0(
  "Supplementary Figure S1. Baseline expression of canonical drug targets, transcriptional assessment of cytotoxicity ",
  "and additional drug-specific analyses. ",
  sprintf("(A) Expression (log2 TPM + 1) of the canonical molecular targets of each substance in media-only control KNS-42 (n = %d) and SVG-A (n = %d) cells; ", nK, nS),
  "points are replicates and bars are means; the dashed line marks TPM = 1. Expressed: + = TPM >= 1 in all replicates; ",
  "(+) = TPM < 1 but >= 10 reads in all replicates; - = not detected (K = KNS-42, S = SVG-A). ",
  "Log2FC and adjusted p compare KNS-42 with SVG-A (DESeq2, Benjamini-Hochberg; * < 0.05, ** < 0.01, *** < 0.001); ",
  "the two lines were sequenced in separate runs, so these comparisons are descriptive. ",
  "The right-hand grid shows regulation of each target by each drug (adjusted p < 0.05, |log2FC| > 1). ",
  "(B) Transcriptional hallmarks of cytotoxicity at 24 h: p53/DNA-damage response, apoptosis effectors, unfolded protein response, ",
  "NRF2-dependent oxidative stress, heat shock, and proliferation/cell cycle genes. Colour shows the median log2 fold-change of each set ",
  "minus that of all other quantified genes; numbers give set genes with adjusted p < 0.05 and |log2FC| > 1; ",
  sprintf("outlined cells differ from background at FDR < 0.05 (two-sided Wilcoxon rank-sum, Benjamini-Hochberg; %d of %d tests). ", n_hit, n_tests),
  "Gene-level values are given in Supplementary Table S2. ",
  "(C) Log2 fold-change of drug-unique genes (columns) in response to every drug (rows), grouped by the drug to which they are unique, ",
  "for KNS-42 (top) and SVG-A (bottom); values capped at +/-5; grey, not quantified. ",
  sprintf("(D) Psilocin-unique genes detected in both cell lines (n = %d): log2 fold-change in KNS-42 against SVG-A, ", nrow(pd_E)),
  "coloured by directional agreement (concordant, green; discordant, red). ",
  "(E) Psilocin-unique metabolite features ranked by log2 fold-change (red, increased; blue, decreased). ",
  "Genes and metabolites: adjusted p < 0.05, |log2FC| > 1. SVG-A morphine was excluded throughout.")
writeLines(strwrap(legend, 100), file.path(out_dir, "FigS1_legend.txt"))
cat("\n", legend, "\n")

if (WRITE_XLSX) write_supplementary(file.path(out_dir, "Supplementary_Tables.xlsx"))

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
message(sprintf("Supplementary Figure S1 written to: %s (%.0f x %.0f mm)",
                out_dir, FIG_W_MM, fig_h * 25.4))
