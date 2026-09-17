

# ---- 1. Settings ----------------------------------------------------------

analysis_dir <- "C:/Users/mrzsjs1/Desktop/Addiction/analysis"
out_dir      <- file.path(analysis_dir, "figures", "Overlap_checks")
PHOS_DIR     <- analysis_dir                        # *_phosphoproteomics_significant.xlsx
METAB_XLSX   <- "C:/Users/mrzsjs1/Desktop/Addiction/analysis/metabolomics final.xlsx"
METAB_SHEET  <- "Filtered"


PTM_TSV         <- "C:/Users/mrzsjs1/Downloads/20241018_112329_P0904_Phospho_Report_PTM.tsv"
PHOS_UNIVERSE_N <- 3992

N_PERM <- 1e6          # permutations per test (sequential sampling makes this fast)
SEED   <- 20260916
PADJ <- 0.05; LFC <- 1; PHOS_P <- 0.05

KNS_DRUGS <- c("ethanol", "THC", "cocaine", "morphine", "psilocin")
SVG_DRUGS <- c("ethanol", "THC", "cocaine", "psilocin")      # SVG-A morphine excluded
tx_file <- function(prefix, drug) {
  if (prefix == "KNS" && drug %in% c("ethanol", "psilocin"))
    return(sprintf("KNS_%s_DESeq2_transcript_results_CORRECTED.csv", drug))
  sprintf("%s_%s_DESeq2_transcript_results.csv", prefix, drug)
}
METAB_CONTRAST <- list(ethanol = c("E","UC"), THC = c("T","EC"), cocaine = c("C","UC"),
                       morphine = c("M","UC"), psilocin = c("P","UC"))


# ---- 2. Packages ----------------------------------------------------------

pkgs <- c("readr", "readxl", "dplyr", "data.table")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing)
suppressPackageStartupMessages({ library(readr); library(readxl); library(dplyr) })
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(SEED)


# ---- 3. Tests ---------------------------------------------------------------

## sequential exact sampler: universes = list of character vectors (detected
## features per drug), sizes = observed list sizes; returns N_PERM intersections
perm_intersection <- function(universes, sizes, n_perm = N_PERM) {
  all_ids <- unique(unlist(universes))
  inU <- lapply(universes, function(u) all_ids %in% u)
  ord <- order(sizes)                                  # smallest list first
  first <- ord[1]
  pool1 <- which(inU[[first]])
  out <- integer(n_perm)
  for (b in seq_len(n_perm)) {
    I <- pool1[sample.int(length(pool1), sizes[first])]
    for (j in ord[-1]) {
      elig <- I[inU[[j]][I]]
      if (!length(elig)) { I <- integer(0); break }
      x <- rhyper(1, length(elig), sum(inU[[j]]) - length(elig), sizes[j])
      I <- if (x == length(elig)) elig else elig[sample.int(length(elig), x)]
    }
    out[b] <- length(I)
  }
  out
}

## exact P(intersection >= obs) for random subsets of a common universe of size N
exact_p <- function(N, sizes, obs) {
  sizes <- sort(sizes)
  d <- numeric(sizes[1] + 1); d[sizes[1] + 1] <- 1
  for (n in sizes[-1]) {
    nd <- numeric(length(d))
    for (m in which(d > 0) - 1) {
      x <- 0:m
      nd[x + 1] <- nd[x + 1] + d[m + 1] * dhyper(x, m, N - m, n)
    }
    d <- nd
  }
  if (obs > length(d) - 1) return(0)
  sum(d[(obs:(length(d) - 1)) + 1])
}

fmt_p <- function(b, n) {
  if (b == 0) sprintf("< %s", formatC(1 / (n + 1), format = "e", digits = 1))
  else formatC((b + 1) / (n + 1), format = "e", digits = 2)
}

run_test <- function(label, sig_sets, universes) {
  ## significant lists restricted to their own detected universe
  sig_sets <- Map(intersect, sig_sets, universes)
  sizes <- lengths(sig_sets)
  shared <- sort(Reduce(intersect, sig_sets))
  obs <- length(shared)
  message(sprintf("%s: sizes %s; observed %d; running %s permutations ...",
                  label, paste(sizes, collapse = "/"), obs, format(N_PERM, big.mark = ",")))
  null <- perm_intersection(universes, sizes)
  b <- sum(null >= obs)
  common <- Reduce(intersect, universes)
  sizes_c <- lengths(lapply(sig_sets, intersect, common))
  data.frame(
    test = label,
    list_sizes = paste(sizes, collapse = "/"),
    universe_sizes = paste(lengths(universes), collapse = "/"),
    observed = obs,
    shared = paste(shared, collapse = ", "),
    perm_expected = round(mean(null), 3),
    perm_null_max = max(null),
    n_perm = N_PERM,
    perm_b = b,
    perm_p = if (b == 0) 1 / (N_PERM + 1) else (b + 1) / (N_PERM + 1),
    perm_p_text = fmt_p(b, N_PERM),
    common_universe = length(common),
    analytic_expected = round(length(common) * prod(sizes_c / length(common)), 3),
    analytic_p = exact_p(length(common), sizes_c, obs),   # shared features lie in every universe
    stringsAsFactors = FALSE)
}


# ---- 4. Transcriptomics ------------------------------------------------------

gene_sets <- function(prefix, drug) {
  d <- read_csv(file.path(analysis_dir, tx_file(prefix, drug)), show_col_types = FALSE,
                col_select = c("log2FoldChange", "padj", "gene_symbol")) %>%
    mutate(gene = toupper(trimws(gene_symbol))) %>%
    filter(!is.na(gene), gene != "", !is.na(padj), !is.na(log2FoldChange))
  list(universe = unique(d$gene),
       sig = unique(d$gene[d$padj < PADJ & abs(d$log2FoldChange) > LFC]))
}
svg <- lapply(setNames(SVG_DRUGS, SVG_DRUGS), function(d) gene_sets("SVG", d))
kns <- lapply(setNames(KNS_DRUGS, KNS_DRUGS), function(d) gene_sets("KNS", d))

res <- list(
  run_test("SVG-A DEGs shared by all four drugs",
           lapply(svg, `[[`, "sig"), lapply(svg, `[[`, "universe")),
  run_test("KNS-42 DEGs shared by all five drugs",
           lapply(kns, `[[`, "sig"), lapply(kns, `[[`, "universe")))


# ---- 5. Phosphoproteomics ----------------------------------------------------

phos_sig <- lapply(setNames(KNS_DRUGS, KNS_DRUGS), function(d) {
  x <- read_excel(file.path(PHOS_DIR, sprintf("%s_phosphoproteomics_significant.xlsx", tolower(d))))
  g <- toupper(trimws(sub("[;,].*$", "", as.character(x$gene[x$p_value < PHOS_P]))))
  unique(g[!is.na(g) & g != ""])
})
phos_universe <- NULL
if (file.exists(PTM_TSV)) {
  hdr <- names(read.delim(PTM_TSV, nrows = 1, check.names = FALSE))
  gcol <- intersect(c("PG.Genes", "PTM.Genes", "EG.Genes", "PG.GeneNames"), hdr)[1]
  if (!is.na(gcol)) {
    g <- data.table::fread(PTM_TSV, select = gcol, sep = "\t", showProgress = FALSE)[[1]]
    g <- toupper(trimws(sub("[;,].*$", "", g)))
    phos_universe <- unique(c(g[!is.na(g) & g != ""], unlist(phos_sig)))
    message("Phospho universe from PTM report (", gcol, "): ", length(phos_universe), " proteins")
  }
}
if (is.null(phos_universe)) {
  ## no gene names available: use PHOS_UNIVERSE_N placeholder identities that
  ## contain every regulated protein
  reg <- unique(unlist(phos_sig))
  phos_universe <- c(reg, sprintf("unregulated_%d", seq_len(PHOS_UNIVERSE_N - length(reg))))
  message("Phospho universe: PHOS_UNIVERSE_N = ", PHOS_UNIVERSE_N, " (PTM report not found)")
}
res[[length(res) + 1]] <- run_test("KNS-42 phosphoproteins shared by all five drugs",
                                   phos_sig, rep(list(phos_universe), length(phos_sig)))


# ---- 6. Metabolomics ---------------------------------------------------------

MF <- read_excel(METAB_XLSX, sheet = METAB_SHEET)
met <- lapply(METAB_CONTRAST, function(cc) {
  d <- data.frame(name = trimws(as.character(MF$Name)),
                  lfc  = -as.numeric(MF[[sprintf("Log2 Fold Change: (%s) / (%s)", cc[2], cc[1])]]),
                  padj =  as.numeric(MF[[sprintf("Adj. P-value: (%s) / (%s)", cc[2], cc[1])]]))
  d <- d[!is.na(d$name) & d$name != "" & !is.na(d$lfc) & !is.na(d$padj), ]
  list(universe = unique(d$name), sig = unique(d$name[d$padj < PADJ & abs(d$lfc) > LFC]))
})
res[[length(res) + 1]] <- run_test("KNS-42 metabolites (named) shared by all five drugs",
                                   lapply(met, `[[`, "sig"), lapply(met, `[[`, "universe"))

## the calculation behind the current text: feature rows as the universe and
## feature counts as list sizes, although the 13 shared entries are names
named <- !is.na(MF$Name) & trimws(as.character(MF$Name)) != ""
feat_sizes <- vapply(METAB_CONTRAST, function(cc) {
  l <- -as.numeric(MF[[sprintf("Log2 Fold Change: (%s) / (%s)", cc[2], cc[1])]])
  p <-  as.numeric(MF[[sprintf("Adj. P-value: (%s) / (%s)", cc[2], cc[1])]])
  sum(named & p < PADJ & abs(l) > LFC, na.rm = TRUE)       # 80/51/72/79/136 in the text
}, numeric(1))
text_calc <- data.frame(
  test = "Metabolites as calculated for the current text (named feature rows vs all feature rows; not consistent with name-level sharing)",
  list_sizes = paste(feat_sizes, collapse = "/"), universe_sizes = as.character(nrow(MF)),
  observed = res[[length(res)]]$observed, shared = "", perm_expected = NA, perm_null_max = NA,
  n_perm = NA, perm_b = NA, perm_p = NA, perm_p_text = "",
  common_universe = nrow(MF),
  analytic_expected = round(nrow(MF) * prod(feat_sizes / nrow(MF)), 3),
  analytic_p = exact_p(nrow(MF), feat_sizes, res[[length(res)]]$observed))


# ---- 7. Report ---------------------------------------------------------------

out <- bind_rows(c(res, list(text_calc)))
write_csv(out, file.path(out_dir, "overlap_permutation_checks.csv"))

cat("\n================ Overlap checks ================\n")
for (i in seq_len(nrow(out))) {
  r <- out[i, ]
  cat(sprintf("\n%s\n  lists %s; detected %s\n  observed %d%s\n",
              r$test, r$list_sizes, r$universe_sizes, r$observed,
              if (nzchar(r$shared)) paste0(" (", r$shared, ")") else ""))
  if (!is.na(r$n_perm))
    cat(sprintf("  permutation: expected %.2f, max in null %d, %d of %s >= observed, p %s\n",
                r$perm_expected, r$perm_null_max, r$perm_b, format(r$n_perm, big.mark = ","),
                r$perm_p_text))
  cat(sprintf("  analytic (common universe %d): expected %.3g, p = %s\n",
              r$common_universe, r$analytic_expected, formatC(r$analytic_p, format = "e", digits = 1)))
}
cat("\nSuggested wording: 'permutation p = X' (or 'p < 1/(N+1)' when no permutation\n",
    "reached the observed overlap), with the number of permutations stated once in Methods.\n", sep = "")
message("\nWritten: ", file.path(out_dir, "overlap_permutation_checks.csv"))
