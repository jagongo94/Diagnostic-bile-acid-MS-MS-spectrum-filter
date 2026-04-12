# Extracting alkamine bile acids MS/MS spectra with diagnostic peaks and the cutoff ratio 
# Extracting alkamine bile acids MS/MS spectra with diagnostic peaks and the cutoff ratio 
# 1. Set working directory
setwd("C:/Users/jagongo/Desktop/Reductive amination project_July 21_2025/Bile acid alkylamine library")

# MS/MS spectra for Lithocholic acid-alkamine conjugates bile acids (Diagnostic peak of 359.2930, 341.2831, 323.2722)
# MS/MS spectra for Lithocholic acid-alkamine bile acids (Diagnostic peak of 359.2930, 341.2831, 323.2722)

# 2. Function to extract full blocks (BEGIN IONS … END IONS) for matching NAME
extract_blocks <- function(file_path, keyword) {
  lines <- readLines(file_path, warn = FALSE)
  
  begin_idx <- grep("^BEGIN IONS", lines)
  end_idx   <- grep("^END IONS", lines)
  
  blocks <- list()
  
  for (i in seq_along(begin_idx)) {
    block_lines <- lines[begin_idx[i]:end_idx[i]]
    name_line <- grep("^NAME=", block_lines, value = TRUE)
    
    if (length(name_line) > 0 && grepl(keyword, name_line, ignore.case = TRUE)) {
      blocks[[length(blocks) + 1]] <- block_lines
    }
  }
  return(blocks)
}

# 3. Function to extract mz/intensity values from a block
extract_peaks <- function(block) {
  peak_lines <- grep("^[0-9]", block, value = TRUE)
  if (length(peak_lines) == 0) return(data.frame(mz=numeric(), intensity=numeric()))
  
  peaks <- do.call(rbind, strsplit(peak_lines, "\\s+"))
  peaks <- peaks[, colSums(!is.na(peaks)) > 0]
  peaks <- data.frame(mz = as.numeric(peaks[,1]), intensity = as.numeric(peaks[,2]))
  return(peaks)
}

# 4. Extract blocks for "Dehydrolithocholicacid"
mgf_file <- "GNPS-ALKYLAMINES-BILE-ACIDS-LIBRARY.mgf"
dhca_blocks <- extract_blocks(mgf_file, "Dehydrolithocholicacid")

# 5. Define target m/z values and ppm tolerance
target_mzs <- c(359.2930, 341.2831, 323.2722)
ppm_tol <- 20

# 6. Prepare results table (wide format)
results <- data.frame(
  NAME=character(),
  SCANS=character(),
  SPECTRUMID=character(),
  intensity_359=numeric(),
  intensity_341=numeric(),
  intensity_323=numeric(),
  stringsAsFactors=FALSE
)

for (block in dhca_blocks) {
  # Extract NAME
  name_line <- grep("^NAME=", block, value = TRUE)
  name_val <- sub("^NAME=", "", name_line)
  
  # Extract SCANS (if exists)
  scans_line <- grep("^SCANS=", block, value = TRUE)
  scans_val <- if(length(scans_line) > 0) sub("^SCANS=", "", scans_line) else NA
  
  # Extract SPECTRUMID (if exists)
  spectrumid_line <- grep("^SPECTRUMID=", block, value = TRUE)
  spectrumid_val <- if(length(spectrumid_line) > 0) sub("^SPECTRUMID=", "", spectrumid_line) else NA
  
  # Extract peaks
  peaks <- extract_peaks(block)
  
  # Initialize intensities as 0
  intensities <- c(0, 0, 0)
  
  if (nrow(peaks) > 0) {
    for (j in seq_along(target_mzs)) {
      target_mz <- target_mzs[j]
      mz_tol <- target_mz * ppm_tol * 1e-6
      match_idx <- which(abs(peaks$mz - target_mz) <= mz_tol)
      
      if (length(match_idx) > 0) {
        # Take the peak with highest intensity if multiple
        intensities[j] <- max(peaks$intensity[match_idx])
      }
    }
  }
  
  # Add row to results
  results <- rbind(results, data.frame(
    NAME = name_val,
    SCANS = scans_val,
    SPECTRUMID = spectrumid_val,
    intensity_359 = intensities[1],
    intensity_341 = intensities[2],
    intensity_323 = intensities[3],
    stringsAsFactors = FALSE
  ))
}

# 7. Save results to CSV
output_csv <- "Dehydrolithocholicacid_peaks_with_scans_spectrumid.csv"
write.csv(results, file = output_csv, row.names = FALSE)

cat("Matrix-style results with SCANS and SPECTRUMID saved to:", output_csv, "\n")



# MS/MS spectra for deoxycholic acid-alkamine conjugates bile acids (Diagnostic peak of 357.2773, 339.2677, 321.2568)
# MS/MS spectra for deoxycholic acid-alkamine conjugates bile acids (Diagnostic peak of 357.2773, 339.2677, 321.2568)

# 2. Function to extract full blocks (BEGIN IONS … END IONS) for matching NAME
extract_blocks <- function(file_path, keyword) {
  lines <- readLines(file_path, warn = FALSE)
  
  begin_idx <- grep("^BEGIN IONS", lines)
  end_idx   <- grep("^END IONS", lines)
  
  blocks <- list()
  
  for (i in seq_along(begin_idx)) {
    block_lines <- lines[begin_idx[i]:end_idx[i]]
    name_line <- grep("^NAME=", block_lines, value = TRUE)
    
    if (length(name_line) > 0 && grepl(keyword, name_line, ignore.case = TRUE)) {
      blocks[[length(blocks) + 1]] <- block_lines
    }
  }
  return(blocks)
}

# 3. Function to extract mz/intensity values from a block
extract_peaks <- function(block) {
  peak_lines <- grep("^[0-9]", block, value = TRUE)
  if (length(peak_lines) == 0) return(data.frame(mz=numeric(), intensity=numeric()))
  
  peaks <- do.call(rbind, strsplit(peak_lines, "\\s+"))
  peaks <- peaks[, colSums(!is.na(peaks)) > 0]
  peaks <- data.frame(mz = as.numeric(peaks[,1]), intensity = as.numeric(peaks[,2]))
  return(peaks)
}

# 4. Extract blocks for "3-keto-deoxycholic acid"
mgf_file <- "GNPS-ALKYLAMINES-BILE-ACIDS-LIBRARY.mgf"
blocks_3keto <- extract_blocks(mgf_file, "3-keto-deoxycholic acid")

# 5. Define target m/z values and ppm tolerance
target_mzs <- c(357.2773, 339.2677, 321.2568)
ppm_tol <- 20

# 6. Prepare results table (wide format)
results <- data.frame(
  NAME=character(),
  SCANS=character(),
  SPECTRUMID=character(),
  intensity_357=numeric(),
  intensity_339=numeric(),
  intensity_321=numeric(),
  stringsAsFactors=FALSE
)

for (block in blocks_3keto) {
  # Extract NAME
  name_line <- grep("^NAME=", block, value = TRUE)
  name_val <- sub("^NAME=", "", name_line)
  
  # Extract SCANS (if exists)
  scans_line <- grep("^SCANS=", block, value = TRUE)
  scans_val <- if(length(scans_line) > 0) sub("^SCANS=", "", scans_line) else NA
  
  # Extract SPECTRUMID (if exists)
  spectrumid_line <- grep("^SPECTRUMID=", block, value = TRUE)
  spectrumid_val <- if(length(spectrumid_line) > 0) sub("^SPECTRUMID=", "", spectrumid_line) else NA
  
  # Extract peaks
  peaks <- extract_peaks(block)
  
  # Initialize intensities as 0
  intensities <- c(0, 0, 0)
  
  if (nrow(peaks) > 0) {
    for (j in seq_along(target_mzs)) {
      target_mz <- target_mzs[j]
      mz_tol <- target_mz * ppm_tol * 1e-6
      match_idx <- which(abs(peaks$mz - target_mz) <= mz_tol)
      
      if (length(match_idx) > 0) {
        # Take the peak with highest intensity if multiple
        intensities[j] <- max(peaks$intensity[match_idx])
      }
    }
  }
  
  # Add row to results
  results <- rbind(results, data.frame(
    NAME = name_val,
    SCANS = scans_val,
    SPECTRUMID = spectrumid_val,
    intensity_357 = intensities[1],
    intensity_339 = intensities[2],
    intensity_321 = intensities[3],
    stringsAsFactors = FALSE
  ))
}

# 7. Save results to CSV
output_csv <- "3keto_deoxycholic_acid_peaks.csv"
write.csv(results, file = output_csv, row.names = FALSE)

cat("Matrix-style results with SCANS and SPECTRUMID saved to:", output_csv, "\n")



# MS/MS spectra for cholic acid-alkamine conjugates bile acids (Diagnostic peak of 355.2626, 337.2531, 319.2416)

# 2. Function to extract full blocks (BEGIN IONS … END IONS) for matching NAME
extract_blocks <- function(file_path, keyword) {
  lines <- readLines(file_path, warn = FALSE)
  
  begin_idx <- grep("^BEGIN IONS", lines)
  end_idx   <- grep("^END IONS", lines)
  
  blocks <- list()
  
  for (i in seq_along(begin_idx)) {
    block_lines <- lines[begin_idx[i]:end_idx[i]]
    name_line <- grep("^NAME=", block_lines, value = TRUE)
    
    if (length(name_line) > 0 && grepl(keyword, name_line, ignore.case = TRUE)) {
      blocks[[length(blocks) + 1]] <- block_lines
    }
  }
  return(blocks)
}

# 3. Function to extract mz/intensity values from a block
extract_peaks <- function(block) {
  peak_lines <- grep("^[0-9]", block, value = TRUE)
  if (length(peak_lines) == 0) return(data.frame(mz=numeric(), intensity=numeric()))
  
  peaks <- do.call(rbind, strsplit(peak_lines, "\\s+"))
  peaks <- peaks[, colSums(!is.na(peaks)) > 0]
  peaks <- data.frame(mz = as.numeric(peaks[,1]), intensity = as.numeric(peaks[,2]))
  return(peaks)
}

# 4. Extract blocks for "3-keto-cholic acid"
mgf_file <- "GNPS-ALKYLAMINES-BILE-ACIDS-LIBRARY.mgf"
blocks_3keto_cholic <- extract_blocks(mgf_file, "3-keto-cholic acid")

# 5. Define new target m/z values and ppm tolerance
target_mzs <- c(355.2626, 337.2531, 319.2416)
ppm_tol <- 20

# 6. Prepare results table (wide format)
results <- data.frame(
  NAME=character(),
  SCANS=character(),
  SPECTRUMID=character(),
  intensity_355=numeric(),
  intensity_337=numeric(),
  intensity_319=numeric(),
  stringsAsFactors=FALSE
)

for (block in blocks_3keto_cholic) {
  # Extract NAME
  name_line <- grep("^NAME=", block, value = TRUE)
  name_val <- sub("^NAME=", "", name_line)
  
  # Extract SCANS (if exists)
  scans_line <- grep("^SCANS=", block, value = TRUE)
  scans_val <- if(length(scans_line) > 0) sub("^SCANS=", "", scans_line) else NA
  
  # Extract SPECTRUMID (if exists)
  spectrumid_line <- grep("^SPECTRUMID=", block, value = TRUE)
  spectrumid_val <- if(length(spectrumid_line) > 0) sub("^SPECTRUMID=", "", spectrumid_line) else NA
  
  # Extract peaks
  peaks <- extract_peaks(block)
  
  # Initialize intensities as 0
  intensities <- c(0, 0, 0)
  
  if (nrow(peaks) > 0) {
    for (j in seq_along(target_mzs)) {
      target_mz <- target_mzs[j]
      mz_tol <- target_mz * ppm_tol * 1e-6
      match_idx <- which(abs(peaks$mz - target_mz) <= mz_tol)
      
      if (length(match_idx) > 0) {
        # Take the peak with highest intensity if multiple
        intensities[j] <- max(peaks$intensity[match_idx])
      }
    }
  }
  
  # Add row to results
  results <- rbind(results, data.frame(
    NAME = name_val,
    SCANS = scans_val,
    SPECTRUMID = spectrumid_val,
    intensity_355 = intensities[1],
    intensity_337 = intensities[2],
    intensity_319 = intensities[3],
    stringsAsFactors = FALSE
  ))
}

# 7. Save results to CSV
output_csv <- "3keto_cholic_acid_peaks.csv"
write.csv(results, file = output_csv, row.names = FALSE)

cat("Matrix-style results with SCANS and SPECTRUMID saved to:", output_csv, "\n")



