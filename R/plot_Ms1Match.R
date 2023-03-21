#' Plots the Biomolecule Isotope Profile on top of the experimental sequence
#'
#' @description Returns a static plot with identified calculated peaks
#'     plotted over the experimental spectrum.
#'
#' @param PeakData (peak_data object) A pspecterlib peak_data object or data.table with "M/Z" and "Intensity". Required.
#' @param Ms1Match (IsoMatchMS1_MatchedPeaks) A IsoMatchMS1_MatchedPeaks object from match_biomolecule_to_ms1. Required.
#' @param ID (numeric) The ID in the IsoMatchMS_MatchedPeaks object to plot. Required.
#' @param Window (numeric) The -/+ m/z value on either side of the matched spectra plot. Default is 5 m/z.
#'
#' @returns (ggplot object) A ggplot object
#'
#' @examples
#' \dontrun{
#'
#' # Run two examples with two charge states
#' MolForms_Test <- calculate_molform(
#'   Biomolecule = c("M.SS[Methyl]S.V", "M.SS[6]S[7].V"),
#'   BioType = "ProForma",
#'   Identifiers = c("Test1", "Test2"),
#'   Charge = 1:2
#' )
#'
#' # Generate some experimental peak data to match
#' PeakData <- pspecterlib::make_peak_data(
#'   MZ = c(147.5684, 148.0699, 148.5708, 149.0721, 149.5731,
#'          294.1296, 295.1325, 296.1343, 297.1369, 298.1390),
#'   Intensity = c(868.3680036, 110.9431876, 18.7179196, 1.7871629, 0.1701294,
#'                 868.3680036, 110.9431876, 18.7179196, 1.7871629, 0.1701294)
#' )
#'
#' # Run algorithm
#' IsoMatch <- match_biomolecule_to_ms1(
#'   PeakData = PeakData,
#'   MolecularFormula = MolForms_Test,
#'   IsotopeMinimum = 3
#' )
#'
#' # Make plot
#' plot_Ms1Match(
#'   PeakData = PeakData,
#'   Ms1Match = IsoMatch,
#'   ID = 1
#' )
#'
#' }
#'
#' @export
plot_Ms1Match <- function(PeakData,
                          Ms1Match,
                          ID,
                          Window = 2) {

  ##################
  ## CHECK INPUTS ##
  ##################

  # Check that peak data is of the appropriate class
  if ("peak_data" %in% class(PeakData) == FALSE) {
    stop("PeakData must be a pspecterlib peak_data object.")
  }

  # Ms1Match should be a IsoMatchMS_MatchedPeaks object
  if ("IsoMatchMS_MatchedPeaks" %in% class(Ms1Match) == FALSE) {
    stop("Ms1Match must be a IsoMatchMS1_MatchedPeaks object from match_biomolecule_to_ms1")
  }

  # ID should be length one
  if (length(ID) != 1) {
    stop("ID should be of length 1.")
  }

  # ID must be an acceptable option
  if (ID %in% Ms1Match$ID == FALSE) {
    stop(paste(ID, "is not a recognized ID."))
  }

  # Window must be numeric
  if (!is.numeric(Window) | length(Window) != 1) {
    stop("Window must be a single numeric.")
  }
  Window <- abs(Window)

  ###################################
  ## MAKE A DATAFRAME FOR PLOTTING ##
  ###################################

  # Change Ms1Match so that it can be used with dplyr functions
  class(PeakData) <- c("data.table", "data.frame")
  class(Ms1Match) <- c("data.table", "data.frame")

  # Filter Ms1 Match to the correct subset
  IDSelection <- ID
  Ms1MatchSub <- Ms1Match %>% dplyr::filter(ID == IDSelection)

  # Adjust PeakData to be within range, rename intensity to abundance
  AdjPeakData <- PeakData %>%
    dplyr::filter(`M/Z` >= min(Ms1MatchSub$`M/Z`) - Window & `M/Z` <= max(Ms1MatchSub$`M/Z`) + Window)


  ###############
  ## MAKE PLOT ##
  ###############

  # Calculate a scaled abundance
  Scale <- Ms1MatchSub[Ms1MatchSub$Abundance == max(Ms1MatchSub$Abundance), "Abundance Experimental"] / max(Ms1MatchSub$Abundance)
  Ms1MatchSub$`Abundance Scaled` <- round(Ms1MatchSub$Abundance * Scale, 4)
  Ms1MatchSub$Matched <- factor(ifelse(!is.na(Ms1MatchSub$`M/Z Experimental`), "Yes", "No"), levels = c("Yes", "No"))

  # Zero fill MS1 match
  MS1 <- data.table::data.table(
    `M/Z` = c(AdjPeakData$`M/Z` - 1e-9, AdjPeakData$`M/Z`, AdjPeakData$`M/Z` + 1e-9),
    Abundance = c(rep(0, nrow(AdjPeakData)), AdjPeakData$Abundance, rep(0, nrow(AdjPeakData)))
  ) %>%
    dplyr::arrange(`M/Z`)

  # Generate the plot
  plot <- ggplot2::ggplot() +
    ggplot2::geom_line(data = MS1, ggplot2::aes(x = `M/Z`, y = Abundance), color = "black", alpha = 0.25) +
    ggplot2::ylim(c(0, max(MS1$Abundance + 1))) +
    ggplot2::geom_point(data = Ms1MatchSub, ggplot2::aes(x = `M/Z`, y = `Abundance Scaled`, color = Matched), size = 3) +
    ggplot2::scale_color_manual(values = c("purple", "orange")) +
    ggplot2::theme_bw() + 
    ggplot2::geom_vline(xintercept = Ms1MatchSub$`Monoisotopic Mass`[1], linetype = "dotted", color = "steelblue", size = 1.5)

  return(plot %>% plotly::ggplotly())

}
