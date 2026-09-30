# Allele frequency calculation and visualization for diploid genotype data.

allele_frequencies <- function(genotypes, marker_names = colnames(genotypes)) {
  if (is.data.frame(genotypes)) {
    genotypes <- as.matrix(genotypes)
  }
  if (!is.matrix(genotypes) || !is.numeric(genotypes)) {
    stop("genotypes must be a numeric matrix or data frame")
  }
  if (nrow(genotypes) == 0L || ncol(genotypes) == 0L) {
    stop("genotypes must contain at least one sample and one marker")
  }
  if (any(genotypes < 0 | genotypes > 2, na.rm = TRUE) ||
      any(genotypes != floor(genotypes), na.rm = TRUE)) {
    stop("genotypes must be coded as 0, 1, or 2, with NA for missing calls")
  }

  if (is.null(marker_names)) {
    marker_names <- paste0("Marker_", seq_len(ncol(genotypes)))
  }
  if (length(marker_names) != ncol(genotypes)) {
    stop("marker_names must have one value per marker")
  }

  called_alleles <- colSums(!is.na(genotypes)) * 2
  alternate_alleles <- colSums(genotypes, na.rm = TRUE)
  data.frame(
    marker = marker_names,
    alternate_allele_count = as.numeric(alternate_alleles),
    called_allele_count = as.numeric(called_alleles),
    frequency = ifelse(called_alleles == 0, NA_real_,
                       alternate_alleles / called_alleles),
    stringsAsFactors = FALSE
  )
}

plot_allele_frequencies <- function(frequencies, main = "Allele frequency",
                                    ylab = "Alternate allele frequency",
                                    col = "#2C7FB8", las = 2, ...) {
  required <- c("marker", "frequency")
  if (!is.data.frame(frequencies) || !all(required %in% names(frequencies))) {
    stop("frequencies must contain marker and frequency columns")
  }
  if (nrow(frequencies) == 0L) {
    stop("frequencies must contain at least one marker")
  }
  if (any(!is.na(frequencies$frequency) &
          (frequencies$frequency < 0 | frequencies$frequency > 1))) {
    stop("frequency values must be between 0 and 1")
  }

  graphics::barplot(
    frequencies$frequency,
    names.arg = frequencies$marker,
    ylim = c(0, 1),
    main = main,
    ylab = ylab,
    col = col,
    las = las,
    ...
  )
}

visualize_allele_frequencies <- function(genotypes, marker_names = colnames(genotypes),
                                         ...) {
  frequencies <- allele_frequencies(genotypes, marker_names)
  plot_allele_frequencies(frequencies, ...)
  invisible(frequencies)
}
