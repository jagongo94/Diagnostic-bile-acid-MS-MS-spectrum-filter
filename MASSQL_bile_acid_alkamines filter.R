# Extracting alkamine bile acids MS/MS spectra with diagnostic peaks and the cutoff ratio 
# Extracting alkamine bile acids MS/MS spectra with diagnostic peaks and the cutoff ratio 
# 1. Set working directory
setwd("C:/Users/jagongo/Desktop/Reductive amination project_July 21_2025/Bile acid alkylamine library")

# MS/MS spectra for Lithocholic acid-alkamine conjugates bile acids (Diagnostic peak of 359.2930, 341.2831, 323.2722)
# MS/MS spectra for Lithocholic acid-alkamine bile acids (Diagnostic peak of 359.2930, 341.2831, 323.2722)

# 1. Set working directory
setwd("C:/Users/jagongo/Desktop/Reductive amination project_July 21_2025/Bile acid alkylamine library")

# 2. Function to extract full blocks (BEGIN IONS … END IONS) for matching NAME
extract_blocks <- function(file_path, keyword) {
  lines <- readLines(file_path, warn = FALSE)
  
  begin_idx <- grep("^BEGIN IONS", lines)
  end_idx   <- grep("^END IONS", lines)
  
  blocks <- list()
  
  for (i in seq_along(begin_idx)) {
    block_lines <- lines[begin_idx[i]:end_idx[i]]
    name_line   <- grep("^NAME=", block_lines, value = TRUE)
    
    if (length(name_line) > 0 && grepl(keyword, name_line, ignore.case = TRUE)) {
      blocks[[length(blocks) + 1]] <- block_lines
    }
  }
  return(blocks)
}

# 3. Function to extract mz/intensity values from a block
extract_peaks <- function(block) {
  peak_lines <- grep("^[0-9]", block, value = TRUE)
  if (length(peak_lines) == 0) return(data.frame(mz = numeric(), intensity = numeric()))
  
  peaks <- do.call(rbind, strsplit(peak_lines, "\\s+"))
  peaks <- peaks[, colSums(!is.na(peaks)) > 0]
  peaks <- data.frame(mz = as.numeric(peaks[, 1]), intensity = as.numeric(peaks[, 2]))
  return(peaks)
}

# 4. Extract blocks for "Dehydrolithocholicacid"
mgf_file    <- "GNPS-ALKYLAMINES-BILE-ACIDS-LIBRARY.mgf"
dhca_blocks <- extract_blocks(mgf_file, "Dehydrolithocholicacid")

# 5. Define target m/z values and ppm tolerance
target_mzs <- c(359.2930, 341.2831, 323.2722)
ppm_tol    <- 20

# 6. Prepare results table
results <- data.frame(
  NAME               = character(),
  SCANS              = character(),
  SPECTRUMID         = character(),
  intensity_359      = numeric(),
  intensity_341      = numeric(),
  intensity_323      = numeric(),
  ratio_359_341      = numeric(),
  ratio_359_323      = numeric(),
  passes_cutoff      = logical(),
  cutoff_reason      = character(),
  stringsAsFactors   = FALSE
)

# 7. Loop over blocks
for (block in dhca_blocks) {
  
  # --- Extract metadata ---
  name_val <- sub("^NAME=", "", grep("^NAME=", block, value = TRUE))
  
  scans_line     <- grep("^SCANS=", block, value = TRUE)
  scans_val      <- if (length(scans_line) > 0) sub("^SCANS=", "", scans_line) else NA
  
  spectrumid_line <- grep("^SPECTRUMID=", block, value = TRUE)
  spectrumid_val  <- if (length(spectrumid_line) > 0) sub("^SPECTRUMID=", "", spectrumid_line) else NA
  
  # --- Extract peaks ---
  peaks <- extract_peaks(block)
  
  # --- Initialize intensities as NA (not found) ---
  intensities <- rep(NA_real_, 3)   # [1]=359, [2]=341, [3]=323
  
  if (nrow(peaks) > 0) {
    for (j in seq_along(target_mzs)) {
      target_mz <- target_mzs[j]
      mz_tol    <- target_mz * ppm_tol * 1e-6
      match_idx <- which(abs(peaks$mz - target_mz) <= mz_tol)
      
      if (length(match_idx) > 0) {
        intensities[j] <- max(peaks$intensity[match_idx])
      }
    }
  }
  
  i359 <- intensities[1]
  i341 <- intensities[2]
  i323 <- intensities[3]
  
  # --- Compute ratios only when BOTH peaks are present AND each >= 1 ---
  ratio_359_341 <- NA_real_
  ratio_359_323 <- NA_real_
  
  if (!is.na(i359) & !is.na(i341) & i359 >= 1 & i341 >= 1) {
    ratio_359_341 <- i359 / i341
  }
  
  if (!is.na(i359) & !is.na(i323) & i359 >= 1 & i323 >= 1) {
    ratio_359_323 <- i359 / i323
  }
  
  # --- Apply cutoff logic ---
  # PASS if:
  #   - intensity_359 is present AND >= 1, AND
  #   - AT LEAST ONE of (341 or 323) is present AND >= 1, AND
  #   - the ratio for that pair is > 1
  
  pair_341_passes <- !is.na(ratio_359_341) & ratio_359_341 > 1
  pair_323_passes <- !is.na(ratio_359_323) & ratio_359_323 > 1
  
  passes_cutoff <- !is.na(i359) & i359 >= 1 & (pair_341_passes | pair_323_passes)
  
  # --- Descriptive reason for transparency ---
  cutoff_reason <- if (is.na(i359) | i359 < 1) {
    "FAIL: intensity_359 absent or < 1"
  } else if (!pair_341_passes & !pair_323_passes) {
    "FAIL: no valid pair ratio > 1 (peaks absent, < 1, or ratio <= 1)"
  } else {
    passing_pairs <- c(
      if (pair_341_passes) "359/341" else NULL,
      if (pair_323_passes) "359/323" else NULL
    )
    paste("PASS via:", paste(passing_pairs, collapse = " & "))
  }
  
  # --- Append row ---
  results <- rbind(results, data.frame(
    NAME          = name_val,
    SCANS         = scans_val,
    SPECTRUMID    = spectrumid_val,
    intensity_359 = i359,
    intensity_341 = i341,
    intensity_323 = i323,
    ratio_359_341 = ratio_359_341,
    ratio_359_323 = ratio_359_323,
    passes_cutoff = passes_cutoff,
    cutoff_reason = cutoff_reason,
    stringsAsFactors = FALSE
  ))
}

# 8. Save FULL results (all spectra with pass/fail column)
write.csv(results,
          file      = "Dehydrolithocholicacid_all_results.csv",
          row.names = FALSE)

# 9. Save FILTERED results (only spectra that pass all cutoffs)
results_filtered <- results[!is.na(results$passes_cutoff) & results$passes_cutoff == TRUE, ]

write.csv(results_filtered,
          file      = "Dehydrolithocholicacid_filtered_ratio_above1.csv",
          row.names = FALSE)

# 10. Summary
cat("================================================\n")
cat("Total spectra extracted   :", nrow(results), "\n")
cat("Spectra passing cutoff    :", nrow(results_filtered), "\n")
cat("Spectra failing cutoff    :", nrow(results) - nrow(results_filtered), "\n")
cat("------------------------------------------------\n")
cat("Full results    ->  Dehydrolithocholicacid_all_results.csv\n")
cat("Filtered results -> Dehydrolithocholicacid_filtered_ratio_above1.csv\n")
cat("================================================\n")








# MS/MS spectra for deoxycholic acid-alkamine conjugates bile acids (Diagnostic peak of 357.2773, 339.2677, 321.2568)
# MS/MS spectra for deoxycholic acid-alkamine conjugates bile acids (Diagnostic peak of 357.2773, 339.2677, 321.2568)

# 1. Set working directory
setwd("C:/Users/jagongo/Desktop/Reductive amination project_July 21_2025/Bile acid alkylamine library")

# 2. Function to extract full blocks (BEGIN IONS … END IONS) for matching NAME
extract_blocks <- function(file_path, keyword) {
  lines <- readLines(file_path, warn = FALSE)
  
  begin_idx <- grep("^BEGIN IONS", lines)
  end_idx   <- grep("^END IONS", lines)
  
  blocks <- list()
  
  for (i in seq_along(begin_idx)) {
    block_lines <- lines[begin_idx[i]:end_idx[i]]
    name_line   <- grep("^NAME=", block_lines, value = TRUE)
    
    if (length(name_line) > 0 && grepl(keyword, name_line, ignore.case = TRUE)) {
      blocks[[length(blocks) + 1]] <- block_lines
    }
  }
  return(blocks)
}

# 3. Function to extract mz/intensity values from a block
extract_peaks <- function(block) {
  peak_lines <- grep("^[0-9]", block, value = TRUE)
  if (length(peak_lines) == 0) return(data.frame(mz = numeric(), intensity = numeric()))
  
  peaks <- do.call(rbind, strsplit(peak_lines, "\\s+"))
  peaks <- peaks[, colSums(!is.na(peaks)) > 0]
  peaks <- data.frame(mz = as.numeric(peaks[, 1]), intensity = as.numeric(peaks[, 2]))
  return(peaks)
}

# 4. Extract blocks for "3-oxodeoxycholicacid"
mgf_file      <- "GNPS-ALKYLAMINES-BILE-ACIDS-LIBRARY.mgf"
oxodca_blocks <- extract_blocks(mgf_file, "3-oxodeoxycholicacid")

# 5. Define target m/z values and ppm tolerance
# 357.2773 = parent diagnostic peak  (replaces 359.2930)
# 339.2667 = -H2O from 357           (replaces 341.2831)
# 319.2561 = -H2O from 339           (replaces 323.2722)
target_mzs <- c(357.2773, 339.2667, 321.2561)
ppm_tol    <- 20

# 6. Prepare results table
results <- data.frame(
  NAME             = character(),
  SCANS            = character(),
  SPECTRUMID       = character(),
  intensity_357    = numeric(),
  intensity_339    = numeric(),
  intensity_321    = numeric(),
  ratio_357_339    = numeric(),
  ratio_357_321    = numeric(),
  passes_cutoff    = logical(),
  cutoff_reason    = character(),
  stringsAsFactors = FALSE
)

# 7. Loop over blocks
for (block in oxodca_blocks) {
  
  # --- Extract metadata ---
  name_val <- sub("^NAME=", "", grep("^NAME=", block, value = TRUE))
  
  scans_line      <- grep("^SCANS=", block, value = TRUE)
  scans_val       <- if (length(scans_line) > 0) sub("^SCANS=", "", scans_line) else NA
  
  spectrumid_line <- grep("^SPECTRUMID=", block, value = TRUE)
  spectrumid_val  <- if (length(spectrumid_line) > 0) sub("^SPECTRUMID=", "", spectrumid_line) else NA
  
  # --- Extract peaks ---
  peaks <- extract_peaks(block)
  
  # --- Initialize intensities as NA (not found) ---
  intensities <- rep(NA_real_, 3)   # [1]=357, [2]=339, [3]=321
  
  if (nrow(peaks) > 0) {
    for (j in seq_along(target_mzs)) {
      target_mz <- target_mzs[j]
      mz_tol    <- target_mz * ppm_tol * 1e-6
      match_idx <- which(abs(peaks$mz - target_mz) <= mz_tol)
      
      if (length(match_idx) > 0) {
        intensities[j] <- max(peaks$intensity[match_idx])
      }
    }
  }
  
  i357 <- intensities[1]
  i339 <- intensities[2]
  i319 <- intensities[3]
  
  # --- Compute ratios only when BOTH peaks are present AND each >= 1 ---
  ratio_357_339 <- NA_real_
  ratio_357_321 <- NA_real_
  
  if (!is.na(i357) & !is.na(i339) & i357 >= 1 & i339 >= 1) {
    ratio_357_339 <- i357 / i339
  }
  
  if (!is.na(i357) & !is.na(i321) & i357 >= 1 & i321 >= 1) {
    ratio_357_321 <- i357 / i321
  }
  
  # --- Apply cutoff logic ---
  # PASS if:
  #   - intensity_357 is present AND >= 1, AND
  #   - AT LEAST ONE of (339 or 319) is present AND >= 1, AND
  #   - the ratio for that pair is > 1
  
  pair_339_passes <- !is.na(ratio_357_339) & ratio_357_339 > 1
  pair_321_passes <- !is.na(ratio_357_321) & ratio_357_321 > 1
  
  passes_cutoff <- !is.na(i357) & i357 >= 1 & (pair_339_passes | pair_321_passes)
  
  # --- Descriptive reason for transparency ---
  cutoff_reason <- if (is.na(i357) | i357 < 1) {
    "FAIL: intensity_357 absent or < 1"
  } else if (!pair_339_passes & !pair_321_passes) {
    "FAIL: no valid pair ratio > 1 (peaks absent, < 1, or ratio <= 1)"
  } else {
    passing_pairs <- c(
      if (pair_339_passes) "357/339" else NULL,
      if (pair_321_passes) "357/321" else NULL
    )
    paste("PASS via:", paste(passing_pairs, collapse = " & "))
  }
  
  # --- Append row ---
  results <- rbind(results, data.frame(
    NAME          = name_val,
    SCANS         = scans_val,
    SPECTRUMID    = spectrumid_val,
    intensity_357 = i357,
    intensity_339 = i339,
    intensity_321 = i321,
    ratio_357_339 = ratio_357_339,
    ratio_357_321 = ratio_357_321,
    passes_cutoff = passes_cutoff,
    cutoff_reason = cutoff_reason,
    stringsAsFactors = FALSE
  ))
}

# 8. Save FULL results (all spectra with pass/fail column)
write.csv(results,
          file      = "3oxodeoxycholicacid_all_results.csv",
          row.names = FALSE)

# 9. Save FILTERED results (only spectra that pass all cutoffs)
results_filtered <- results[!is.na(results$passes_cutoff) & results$passes_cutoff == TRUE, ]

write.csv(results_filtered,
          file      = "3oxodeoxycholicacid_filtered_ratio_above1.csv",
          row.names = FALSE)

# 10. Summary
cat("================================================\n")
cat("Compound          : 3-Oxodeoxycholic acid\n")
cat("Diagnostic peaks  : 357.2773 | 339.2667 | 321.2561\n")
cat("------------------------------------------------\n")
cat("Total spectra extracted   :", nrow(results), "\n")
cat("Spectra passing cutoff    :", nrow(results_filtered), "\n")
cat("Spectra failing cutoff    :", nrow(results) - nrow(results_filtered), "\n")
cat("------------------------------------------------\n")
cat("Full results     ->  3oxodeoxycholicacid_all_results.csv\n")
cat("Filtered results ->  3oxodeoxycholicacid_filtered_ratio_above1.csv\n")
cat("================================================\n")








# MS/MS spectra for cholic acid-alkamine conjugates bile acids (Diagnostic peak of 355.2626, 337.2531, 319.2416)

# 1. Set working directory
setwd("C:/Users/jagongo/Desktop/Reductive amination project_July 21_2025/Bile acid alkylamine library")

# 2. Function to extract full blocks (BEGIN IONS … END IONS) for matching NAME
extract_blocks <- function(file_path, keyword) {
  lines <- readLines(file_path, warn = FALSE)
  
  begin_idx <- grep("^BEGIN IONS", lines)
  end_idx   <- grep("^END IONS", lines)
  
  blocks <- list()
  
  for (i in seq_along(begin_idx)) {
    block_lines <- lines[begin_idx[i]:end_idx[i]]
    name_line   <- grep("^NAME=", block_lines, value = TRUE)
    
    if (length(name_line) > 0 && grepl(keyword, name_line, ignore.case = TRUE)) {
      blocks[[length(blocks) + 1]] <- block_lines
    }
  }
  return(blocks)
}

# 3. Function to extract mz/intensity values from a block
extract_peaks <- function(block) {
  peak_lines <- grep("^[0-9]", block, value = TRUE)
  if (length(peak_lines) == 0) return(data.frame(mz = numeric(), intensity = numeric()))
  
  peaks <- do.call(rbind, strsplit(peak_lines, "\\s+"))
  peaks <- peaks[, colSums(!is.na(peaks)) > 0]
  peaks <- data.frame(mz = as.numeric(peaks[, 1]), intensity = as.numeric(peaks[, 2]))
  return(peaks)
}

# 4. Extract blocks for "3-oxocholicacid"
mgf_file       <- "GNPS-ALKYLAMINES-BILE-ACIDS-LIBRARY.mgf"
oxoca_blocks   <- extract_blocks(mgf_file, "3-oxocholicacid")

# 5. Define target m/z values and ppm tolerance
# 355.2616 = parent diagnostic peak  (replaces 359.2930)
# 337.2510 = -H2O from 355           (replaces 341.2831)
# 319.2405 = -H2O from 337           (replaces 323.2722)
target_mzs <- c(355.2616, 337.2510, 319.2405)
ppm_tol    <- 20

# 6. Prepare results table
results <- data.frame(
  NAME             = character(),
  SCANS            = character(),
  SPECTRUMID       = character(),
  intensity_355    = numeric(),
  intensity_337    = numeric(),
  intensity_319    = numeric(),
  ratio_355_337    = numeric(),
  ratio_355_319    = numeric(),
  passes_cutoff    = logical(),
  cutoff_reason    = character(),
  stringsAsFactors = FALSE
)

# 7. Loop over blocks
for (block in oxoca_blocks) {
  
  # --- Extract metadata ---
  name_val <- sub("^NAME=", "", grep("^NAME=", block, value = TRUE))
  
  scans_line      <- grep("^SCANS=", block, value = TRUE)
  scans_val       <- if (length(scans_line) > 0) sub("^SCANS=", "", scans_line) else NA
  
  spectrumid_line <- grep("^SPECTRUMID=", block, value = TRUE)
  spectrumid_val  <- if (length(spectrumid_line) > 0) sub("^SPECTRUMID=", "", spectrumid_line) else NA
  
  # --- Extract peaks ---
  peaks <- extract_peaks(block)
  
  # --- Initialize intensities as NA (not found) ---
  intensities <- rep(NA_real_, 3)   # [1]=355, [2]=337, [3]=319
  
  if (nrow(peaks) > 0) {
    for (j in seq_along(target_mzs)) {
      target_mz <- target_mzs[j]
      mz_tol    <- target_mz * ppm_tol * 1e-6
      match_idx <- which(abs(peaks$mz - target_mz) <= mz_tol)
      
      if (length(match_idx) > 0) {
        intensities[j] <- max(peaks$intensity[match_idx])
      }
    }
  }
  
  i355 <- intensities[1]
  i337 <- intensities[2]
  i319 <- intensities[3]
  
  # --- Compute ratios only when BOTH peaks are present AND each >= 1 ---
  ratio_355_337 <- NA_real_
  ratio_355_319 <- NA_real_
  
  if (!is.na(i355) & !is.na(i337) & i355 >= 1 & i337 >= 1) {
    ratio_355_337 <- i355 / i337
  }
  
  if (!is.na(i355) & !is.na(i319) & i355 >= 1 & i319 >= 1) {
    ratio_355_319 <- i355 / i319
  }
  
  # --- Apply cutoff logic ---
  # PASS if:
  #   - intensity_355 is present AND >= 1, AND
  #   - AT LEAST ONE of (337 or 319) is present AND >= 1, AND
  #   - the ratio for that pair is > 1
  
  pair_337_passes <- !is.na(ratio_355_337) & ratio_355_337 > 1
  pair_319_passes <- !is.na(ratio_355_319) & ratio_355_319 > 1
  
  passes_cutoff <- !is.na(i355) & i355 >= 1 & (pair_337_passes | pair_319_passes)
  
  # --- Descriptive reason for transparency ---
  cutoff_reason <- if (is.na(i355) | i355 < 1) {
    "FAIL: intensity_355 absent or < 1"
  } else if (!pair_337_passes & !pair_319_passes) {
    "FAIL: no valid pair ratio > 1 (peaks absent, < 1, or ratio <= 1)"
  } else {
    passing_pairs <- c(
      if (pair_337_passes) "355/337" else NULL,
      if (pair_319_passes) "355/319" else NULL
    )
    paste("PASS via:", paste(passing_pairs, collapse = " & "))
  }
  
  # --- Append row ---
  results <- rbind(results, data.frame(
    NAME          = name_val,
    SCANS         = scans_val,
    SPECTRUMID    = spectrumid_val,
    intensity_355 = i355,
    intensity_337 = i337,
    intensity_319 = i319,
    ratio_355_337 = ratio_355_337,
    ratio_355_319 = ratio_355_319,
    passes_cutoff = passes_cutoff,
    cutoff_reason = cutoff_reason,
    stringsAsFactors = FALSE
  ))
}

# 8. Save FULL results (all spectra with pass/fail column)
write.csv(results,
          file      = "3oxocholicacid_all_results.csv",
          row.names = FALSE)

# 9. Save FILTERED results (only spectra that pass all cutoffs)
results_filtered <- results[!is.na(results$passes_cutoff) & results$passes_cutoff == TRUE, ]

write.csv(results_filtered,
          file      = "3oxocholicacid_filtered_ratio_above1.csv",
          row.names = FALSE)

# 10. Summary
cat("================================================\n")
cat("Compound          : 3-Oxocholic acid\n")
cat("Diagnostic peaks  : 355.2616 | 337.2510 | 319.2405\n")
cat("------------------------------------------------\n")
cat("Total spectra extracted   :", nrow(results), "\n")
cat("Spectra passing cutoff    :", nrow(results_filtered), "\n")
cat("Spectra failing cutoff    :", nrow(results) - nrow(results_filtered), "\n")
cat("------------------------------------------------\n")
cat("Full results     ->  3oxocholicacid_all_results.csv\n")
cat("Filtered results ->  3oxocholicacid_filtered_ratio_above1.csv\n")
cat("================================================\n")
