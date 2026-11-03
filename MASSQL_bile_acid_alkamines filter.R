# =============================================================================
# Diagnostic bile acid MS/MS spectrum filter
# Extracts bile acid–alkylamine conjugate MS/MS spectra from an MGF library and
# keeps those whose diagnostic core ion is present and more intense than at
# least one of its water-loss fragments.
# =============================================================================

# ---- Settings ---------------------------------------------------------------
setwd("C:/Users/jagongo/Desktop/Reductive amination project_July 21_2025/Bile acid alkylamine library")
mgf_file <- "GNPS-ALKYLAMINES-BILE-ACIDS-LIBRARY.mgf"
ppm_tol  <- 20

# One row per bile acid core.
# keyword = regular expression matched (case-insensitive) against NAME=
# ions    = diagnostic ion, -H2O, -2H2O
bile_acids <- list(
  list(label   = "Dehydrolithocholicacid",
       keyword = "dehydrolithocholic ?acid",
       ions    = c(359.2945, 341.2839, 323.2733)),  # C24H39O2+, C24H37O+, C24H35+
  list(label   = "3oxodeoxycholicacid",
       keyword = "3-(keto|oxo)-?deoxycholic ?acid",
       ions    = c(357.2788, 339.2682, 321.2577)),  # C24H37O2+, C24H35O+, C24H33+
  list(label   = "3oxocholicacid",
       keyword = "3-(keto|oxo)-?cholic ?acid",
       ions    = c(355.2632, 337.2526, 319.2420))   # C24H35O2+, C24H33O+, C24H31+
)

# ---- Read the MGF once ------------------------------------------------------
read_mgf_blocks <- function(file_path) {
  lines     <- readLines(file_path, warn = FALSE)
  begin_idx <- grep("^BEGIN IONS", lines)
  end_idx   <- grep("^END IONS", lines)
  lapply(seq_along(begin_idx), function(i) lines[begin_idx[i]:end_idx[i]])
}

get_field <- function(block, key) {
  x <- grep(paste0("^", key, "="), block, value = TRUE)
  if (length(x) > 0) sub(paste0("^", key, "="), "", x[1]) else NA_character_
}

extract_peaks <- function(block) {
  peak_lines <- grep("^[0-9]", block, value = TRUE)
  if (length(peak_lines) == 0) return(data.frame(mz = numeric(), intensity = numeric()))
  parts <- strsplit(trimws(peak_lines), "\\s+")
  data.frame(mz        = as.numeric(vapply(parts, `[`, "", 1)),
             intensity = as.numeric(vapply(parts, `[`, "", 2)))
}

# ---- Filter one bile acid core ----------------------------------------------
filter_bile_acid <- function(all_blocks, label, keyword, ions, ppm_tol) {
  names_all <- vapply(all_blocks, get_field, "", key = "NAME")
  blocks    <- all_blocks[grepl(keyword, names_all, ignore.case = TRUE)]
  tag       <- round(ions)                       # e.g. 359, 341, 323

  rows <- lapply(blocks, function(block) {
    peaks <- extract_peaks(block)
    inten <- vapply(ions, function(t) {
      hit <- which(abs(peaks$mz - t) <= t * ppm_tol * 1e-6)
      if (length(hit) > 0) max(peaks$intensity[hit]) else NA_real_
    }, numeric(1))

    ok    <- !is.na(inten) & inten >= 1
    r1    <- if (ok[1] && ok[2]) inten[1] / inten[2] else NA_real_
    r2    <- if (ok[1] && ok[3]) inten[1] / inten[3] else NA_real_
    p1    <- !is.na(r1) && r1 > 1
    p2    <- !is.na(r2) && r2 > 1
    pass  <- ok[1] && (p1 || p2)

    reason <- if (!ok[1]) {
      sprintf("FAIL: intensity_%d absent or < 1", tag[1])
    } else if (!p1 && !p2) {
      "FAIL: no valid pair ratio > 1 (peaks absent, < 1, or ratio <= 1)"
    } else {
      paste("PASS via:", paste(c(if (p1) sprintf("%d/%d", tag[1], tag[2]),
                                 if (p2) sprintf("%d/%d", tag[1], tag[3])), collapse = " & "))
    }

    out <- data.frame(NAME = get_field(block, "NAME"),
                      SCANS = get_field(block, "SCANS"),
                      SPECTRUMID = get_field(block, "SPECTRUMID"),
                      i1 = inten[1], i2 = inten[2], i3 = inten[3],
                      r1 = r1, r2 = r2,
                      passes_cutoff = pass, cutoff_reason = reason,
                      stringsAsFactors = FALSE)
    names(out)[4:8] <- c(sprintf("intensity_%d", tag),
                         sprintf("ratio_%d_%d", tag[1], tag[2]),
                         sprintf("ratio_%d_%d", tag[1], tag[3]))
    out
  })
  results <- if (length(rows) > 0) do.call(rbind, rows) else data.frame()
  filtered <- if (nrow(results) > 0) results[results$passes_cutoff, ] else results

  all_csv  <- paste0(label, "_all_results.csv")
  filt_csv <- paste0(label, "_filtered_ratio_above1.csv")
  write.csv(results,  all_csv,  row.names = FALSE)
  write.csv(filtered, filt_csv, row.names = FALSE)

  cat("================================================\n")
  cat("Compound          :", label, "\n")
  cat("Diagnostic peaks  :", paste(sprintf("%.4f", ions), collapse = " | "), "\n")
  cat("------------------------------------------------\n")
  cat("Total spectra extracted   :", nrow(results), "\n")
  cat("Spectra passing cutoff    :", nrow(filtered), "\n")
  cat("Spectra failing cutoff    :", nrow(results) - nrow(filtered), "\n")
  cat("Full results     -> ", all_csv, "\n")
  cat("Filtered results -> ", filt_csv, "\n")
  invisible(results)
}

# ---- Run --------------------------------------------------------------------
all_blocks <- read_mgf_blocks(mgf_file)
for (ba in bile_acids) {
  filter_bile_acid(all_blocks, ba$label, ba$keyword, ba$ions, ppm_tol)
}
cat("================================================\n")
