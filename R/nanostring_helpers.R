# =============================================================================
# nanostring_helpers.R
#
# Helper functions for NanoString nCounter normalization in the Fungal AST
# pipeline. Sourced by 01_normalize_nanostring.R.
#
# Modernized from the lab's original `HeatMapFuncs.R`.
#
# -----------------------------------------------------------------------------
# SCOPE OF CHANGES vs. the original HeatMapFuncs.R
#
# The numerical behaviour of every function is IDENTICAL to the original.
# Changes are confined to robustness and readability:
#
#   1. `dplyr::filter(grepl(..., rownames(x)))` replaced with base-R row
#      subsetting. dplyr's row-name handling for base data frames has changed
#      across versions; base subsetting is guaranteed to preserve row names,
#      which the rest of the pipeline depends on absolutely (probe identity is
#      carried entirely in row names). Same rows selected, same order.
#
#   2. `removeCtrlMaxCoV()` now takes `limitCoV` as a REQUIRED argument rather
#      than reading it from the global environment. This is the one signature
#      change in the file: an older script calling `removeCtrlMaxCoV(x)` will
#      now stop with an explanatory error instead of silently using a value it
#      did not choose. 01_normalize_nanostring.R passes the same value the
#      original set globally, so pipeline results are unaffected.
#
#   3. `drop = FALSE` added where subsetting could silently collapse a
#      single-row data frame into a vector.
#
#   4. Documentation and inline comments added throughout.
#
# Function names, argument names, and argument order are UNCHANGED, so any
# existing lab script that sources this file will continue to work.
# =============================================================================


# -----------------------------------------------------------------------------
# Shared colour scheme (kept from the original file)
# -----------------------------------------------------------------------------

#' Blue-white-red diverging palette used for log2 fold-change heatmaps.
myPalette <- colorRampPalette(c("blue", "white", "red"))(n = 299)


# -----------------------------------------------------------------------------
# Small predicates for locating special values in a data frame
#
# These return a logical MATRIX the same shape as the input, so they can be
# used directly on the left-hand side of an assignment:
#     df[is.na.df(df)] <- 0
#
# Idiom adapted from:
# https://stackoverflow.com/questions/18142117/
# -----------------------------------------------------------------------------

#' Logical matrix marking infinite values in a data frame.
#'
#' Infinite values typically appear when a control probe fails and a fold
#' change is computed against a zero denominator.
#'
#' @param x A data frame.
#' @return A logical matrix with the same dimensions as `x`.
is.inf.df <- function(x) do.call(cbind, lapply(x, is.infinite))

#' Logical matrix marking NA values in a data frame.
#'
#' @param x A data frame.
#' @return A logical matrix with the same dimensions as `x`.
is.na.df <- function(x) do.call(cbind, lapply(x, is.na))

#' Logical matrix marking zero values in a data frame.
#'
#' @param x A data frame.
#' @return A logical matrix with the same dimensions as `x`.
is.zero.df <- function(x) do.call(cbind, lapply(x, function(y) y == 0))

#' Logical matrix marking negative values in a data frame.
#'
#' @param x A data frame.
#' @param is_zero_or_less If TRUE, mark values <= 0 rather than < 0.
#' @return A logical matrix with the same dimensions as `x`.
is.neg.df <- function(x, is_zero_or_less = FALSE) {
  if (!is_zero_or_less) {
    do.call(cbind, lapply(x, function(y) y < 0))
  } else {
    do.call(cbind, lapply(x, function(y) y <= 0))
  }
}


# -----------------------------------------------------------------------------
# Internal: row selection by probe-name pattern
#
# Replaces the original `dplyr::filter(grepl(pattern, rownames(myData)))`.
# Always returns a data frame and always preserves row names.
# -----------------------------------------------------------------------------

#' Select rows whose row name matches (or does not match) a regex.
#'
#' @param myData A data frame whose row names are probe identifiers.
#' @param pattern A regular expression matched against `rownames(myData)`.
#' @param invert If TRUE, return the rows that do NOT match.
#' @return A data frame (row names preserved, original order preserved).
#' @keywords internal
.rowsMatching <- function(myData, pattern, invert = FALSE) {
  hit <- grepl(pattern, rownames(myData))
  if (invert) hit <- !hit
  myData[hit, , drop = FALSE]
}


# -----------------------------------------------------------------------------
# Core normalization functions
# -----------------------------------------------------------------------------

#' Geometric mean of every column.
#'
#' The geometric mean is used rather than the arithmetic mean because nCounter
#' counts are multiplicative in nature: a lane that ran "hot" scales every
#' probe by roughly the same factor.
#'
#' @param myData A data frame or matrix of strictly positive counts.
#' @param na.rm If TRUE, NA values are treated as log(1) = 0 before averaging,
#'   i.e. they pull the geometric mean toward 1 rather than propagating NA.
#'   NOTE: this is NOT the same as excluding them from the mean; it is the
#'   original behaviour and is preserved deliberately.
#' @return A named numeric vector, one geometric mean per column.
colGeoMeans <- function(myData, na.rm = FALSE) {
  if (na.rm) {
    myData <- log(myData)
    myData[is.na(myData)] <- 0
    geoMeanData <- exp(colMeans(myData))
  } else {
    geoMeanData <- exp(colMeans(log(myData)))
  }
  return(geoMeanData)
}


#' Positive control correction (lane-to-lane scaling).
#'
#' Every cartridge lane is spiked with the same six ERCC positive controls at
#' known concentrations. Differences in their geometric mean between lanes
#' therefore reflect technical variation (hybridization efficiency, binding
#' density, scan quality) rather than biology. Each lane is divided by its own
#' scale factor so that all lanes share the average positive-control level.
#'
#' Applied to ALL rows, including the controls themselves.
#'
#' @param myData A data frame of counts; row names must contain "POS" for the
#'   positive control probes.
#' @return A data frame the same shape as `myData`, lane-scaled.
posCtrlCorr <- function(myData) {
  # Geometric mean of the POS probes within each lane
  posCtrlGeoMeans <- colGeoMeans(.rowsMatching(myData, "POS"))

  # Scale each lane's geometric mean by the average across lanes.
  # A lane that ran hot gets colScale > 1 and is divided down.
  colScale <- posCtrlGeoMeans / mean(posCtrlGeoMeans)

  posCtrlCorrData <- sweep(myData, 2, colScale, FUN = "/")
  return(posCtrlCorrData)
}


#' Negative control correction (background subtraction).
#'
#' The mean of the negative control probes in a lane estimates that lane's
#' non-specific binding background. That mean is subtracted from every
#' non-control probe.
#'
#' IMPORTANT: the returned data frame EXCLUDES the NEG and POS rows. The
#' calling script re-attaches them with rbind() afterwards. This is the
#' original behaviour.
#'
#' @param myData A data frame of counts; row names must contain "NEG"/"POS"
#'   for the control probes.
#' @param round.to.zero If TRUE, clamp resulting negative values to 0.
#' @return A data frame of background-subtracted non-control probes.
negCtrlCorr <- function(myData, round.to.zero = FALSE) {
  # Per-lane background estimate
  negCtrlMeans <- colMeans(.rowsMatching(myData, "NEG"), na.rm = TRUE)

  # Subtract background from the non-control probes only
  negCtrlCorrData <- sweep(
    .rowsMatching(myData, "NEG|POS", invert = TRUE),
    2, negCtrlMeans
  )

  if (round.to.zero) {
    negCtrlCorrData[is.neg.df(negCtrlCorrData)] <- 0
  }

  return(negCtrlCorrData)
}


#' Negative control correction with floor replacement (alternative method).
#'
#' Variant of `negCtrlCorr()` that floors low values at the background level
#' instead of allowing them to fall below it. With
#' `background.subtraction = TRUE` a probe must exceed roughly twice the
#' background to register any signal at all.
#'
#' NOT used by the current pipeline. Retained from the original file so that
#' older analyses remain reproducible.
#'
#' @param myData A data frame of counts.
#' @param background.subtraction If TRUE, subtract background first, then floor.
#'   If FALSE, only floor.
#' @return A data frame of non-control probes.
negCtrl_wReplacement <- function(myData, background.subtraction = FALSE) {
  negCtrlMeans <- colMeans(.rowsMatching(myData, "NEG"))

  negCtrlCorrData <- .rowsMatching(myData, "NEG|POS", invert = TRUE)

  if (background.subtraction) {
    negCtrlCorrData <- sweep(negCtrlCorrData, 2, negCtrlMeans)
  }

  # Floor each column at its own background value
  for (x in seq_along(negCtrlMeans)) {
    negCtrlCorrData[x][negCtrlCorrData[x] < negCtrlMeans[x]] <- negCtrlMeans[x]
  }

  return(negCtrlCorrData)
}


#' Flag (and optionally drop) probes that fall below a count threshold.
#'
#' A probe is considered failed if ANY sample reads below `minVal`. Failed
#' values are set to NA; with `rm = TRUE` the whole probe row is then dropped,
#' so a probe is only kept if it is reliable in EVERY sample.
#'
#' @param myData A data frame of counts.
#' @param minVal Minimum acceptable count. Typically ~100 for control probes
#'   and ~40 for response probes, though this pipeline uses lower values.
#' @param rm If TRUE, drop failing probe rows entirely (via na.omit).
#' @param firstRow Apply the threshold to the FIRST COLUMN only. Note the
#'   argument name says "row" but the original code indexes `myData[1]`, which
#'   in base R selects the first COLUMN. Name kept for backward compatibility.
#'   In this pipeline the first column is the untreated sample, so this option
#'   means "only require the untreated condition to pass".
#' @param na.rm If TRUE, protect pre-existing NA values from removal by
#'   temporarily converting them to Inf, then restoring them to NA on exit.
#' @return A data frame, with failing values set to NA and (if `rm`) failing
#'   rows removed.
removeFailedProbes <- function(myData, minVal = 50, rm = TRUE,
                               firstRow = FALSE, na.rm = FALSE) {

  # Shield values that were already NA so na.omit() does not drop them for a
  # reason unrelated to the count threshold.
  if (na.rm) {
    myData[is.na.df(myData)] <- Inf
  }

  if (firstRow) {
    # First column only (see @param note above)
    myData[1][(myData[1]) < minVal] <- NA
  } else {
    myData[myData < minVal] <- NA
  }

  if (rm) {
    myData <- na.omit(myData)
    if (na.rm) {
      myData[is.inf.df(myData)] <- NA
    }
    return(myData)
  }

  if (na.rm) {
    myData[is.inf.df(myData)] <- NA
  }
  return(myData)
}


#' Coefficient of variation for each control probe across samples.
#'
#' Each control probe is first expressed as a ratio to the geometric mean of
#' all control probes in its own lane. This removes lane-level scaling, so the
#' remaining spread reflects how inconsistently that particular probe behaves.
#' A high CoV means the probe is a poor normalizer.
#'
#' @param ctrlProbes A data frame of control probe counts (probes x samples).
#' @return A named numeric vector of CoVs, one per probe.
calcCtrlCoVs <- function(ctrlProbes) {
  ctrlProbeRatio      <- sweep(ctrlProbes, 2, colGeoMeans(ctrlProbes), FUN = "/")
  ctrlProbeRatioSDs   <- apply(ctrlProbeRatio, 1, sd)
  ctrlProbeRatioMeans <- rowMeans(ctrlProbeRatio)
  ctrlProbeRatioCoVs  <- ctrlProbeRatioSDs / ctrlProbeRatioMeans
  return(ctrlProbeRatioCoVs)
}


#' One iteration of greedy control-probe pruning.
#'
#' Computes the CoV of every control probe and, if the worst one exceeds
#' `limitCoV`, removes it. Called repeatedly by the driver script until the
#' worst remaining CoV is within tolerance.
#'
#' Note on the return value: element 3 holds the CoVs of the probe set as it
#' was BEFORE this call's removal. The driver script uses `max()` of that
#' vector as the loop condition, which means the loop performs one final
#' pass that removes nothing and then exits. This is correct and intended.
#'
#' @param ctrlProbesCurrent Data frame of control probes still in the running.
#' @param probeToRemove Sentinel default; leave as "none".
#' @param limitCoV Maximum tolerated coefficient of variation. REQUIRED.
#'   The original read this from the global environment. It is now a mandatory
#'   argument with no default, deliberately: a default would silently override
#'   a caller that had set a different global, which is worse than an error.
#'   Older scripts calling `removeCtrlMaxCoV(x)` must be updated to
#'   `removeCtrlMaxCoV(x, limitCoV = limitCoV)`.
#' @return A list of three elements:
#'   \enumerate{
#'     \item the updated control probe data frame,
#'     \item the name of the probe removed (or "none"),
#'     \item the vector of CoVs computed this iteration.
#'   }
removeCtrlMaxCoV <- function(ctrlProbesCurrent, probeToRemove = "none",
                             limitCoV) {

  if (missing(limitCoV)) {
    stop("removeCtrlMaxCoV() now requires an explicit `limitCoV` argument. ",
         "It used to be read from the global environment. ",
         "Call it as removeCtrlMaxCoV(x, limitCoV = limitCoV).")
  }

  # CoV of each probe, measured as its ratio to the lane geometric mean
  ctrlProbeRatioCoVs <- calcCtrlCoVs(ctrlProbesCurrent)

  maxCoV_iter <- max(ctrlProbeRatioCoVs)

  if (maxCoV_iter > limitCoV) {
    probeToRemove <- names(which.max(ctrlProbeRatioCoVs))
  }

  # When probeToRemove is still "none" this removes nothing, because no probe
  # is literally named "none".
  ctrlProbesNext <- ctrlProbesCurrent[
    rownames(ctrlProbesCurrent) != probeToRemove, , drop = FALSE
  ]

  to_return <- list(ctrlProbesNext, probeToRemove, ctrlProbeRatioCoVs)
  return(to_return)
}


#' Per-lane detection floor: six standard deviations of the negative controls.
#'
#' Values below this are indistinguishable from background noise, so the
#' driver script clamps them up to this floor rather than letting near-zero
#' counts produce wild fold changes.
#'
#' @param rdata A data frame of counts; row names must contain "NEG_" for the
#'   negative control probes.
#' @return A named numeric vector, one floor value per lane/column.
minValue6xStDev <- function(rdata) {
  data  <- .rowsMatching(rdata, "NEG_")
  stDev <- apply(data, 2, sd, na.rm = TRUE) * 6
  return(stDev)
}


#' Replace near-zero values in a vector with its first element.
#'
#' Substitutes the untreated value (assumed to be element 1) into any treated
#' position reading below 2, so the probe registers no change rather than
#' being discarded.
#'
#' NOT used by the current pipeline. Retained from the original file.
#'
#' @param vctr A numeric vector whose first element is the untreated value.
#' @return The vector with low treated values replaced.
subLowValues <- function(vctr) {
  len <- length(vctr)
  vctr[2:len][vctr[2:len] < 2] <- vctr[1]
  vctr <- unlist(vctr)
  return(vctr)
}
