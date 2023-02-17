# Get downloads folder
#' @export
.getDownloadsFolder <- function(.test_mode = FALSE) {
  if (Sys.info()['sysname'] == "Windows" | .test_mode) {
    folder <- dirname("~")
    return(folder)
  } else {
    folder <- path.expand("~")
    folder <- file.path(folder, "Downloads")
    folder <- paste0(folder, .Platform$file.sep)
    return(folder)
  }
}

#' Generate a trelliscope display for all Proteoform to MS1 Matches
#'
#' @description [Trelliscope](https://github.com/hafen/trelliscopejs) allows
#'    users to filter and sort plots based on cognostic variables. For this dataset,
#'    the plots are generated with the plot_Ms1Match function, and the cognostics
#'    are: Protein, Absolute Relative Error, Correlation, Charge, Proteoform, and ID.
#'
#' @param PeakData A pspecterlib peak_data object or data.table with "M/Z" and "Intensity". Required.
#' @param Ms1Match A ProteoMatch_MatchedPeaks class object from match_full_seq_ms1. Required.
#' @param Path The base directory of the trelliscope application. Default is Downloads/Ms1Match.
#' @param MinCorrelationScore The minimum correlation score to plot. Default is 0.7.
#' @param Window The -/+ m/z value on either side of the matched spectra plot. Default is 2 m/z.
#'
#' @returns An html trelliscope display
#'
#' @examples
#' \dontrun{
#'
#' # Run two examples with two charge states
#' MolForms_Test <- calculate_molform(
#'   Proteoform = c("M.SS[Methyl]S.V", "M.SS[6]S[7].V"),
#'   Protein = c("Test1", "Test2"),
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
#' ProteoMatch <- match_proteoform_to_ms1(
#'   PeakData = PeakData,
#'   MolecularFormula = MolForms_Test,
#'   IsotopeRange = c(3, 20)
#' )
#'
#' # Make the trelliscope display
#' proteomatch_trelliscope(PeakData = PeakData, Ms1Match = ProteoMatch)
#'
#' }
#'
#'
#' @export
proteomatch_trelliscope <- function(PeakData,
                                    Ms1Match,
                                    Path = file.path(.getDownloadsFolder(), "Ms1Match", "Trelliscope"),
                                    MinCorrelationScore = 0.7,
                                    Window = 2) {

  ##################
  ## CHECK INPUTS ##
  ##################

  # Check that peak data is of the appropriate class
  if ("peak_data" %in% class(PeakData) == FALSE) {
    stop("PeakData must be a pspecterlib peak_data object.")
  }

  # Ms1Match should be a ProteoMatch_MatchedPeaks object
  if ("ProteoMatch_MatchedPeaks" %in% class(Ms1Match) == FALSE) {
    stop("Ms1Match must be a ProteoMatch_MatchedPeaks object from match_proteoform_to_ms1")
  }

  # Check that minimum correlation score is within 0-1
  if (!is.numeric(MinCorrelationScore)) {
    stop("MinCorrelationScore must be a numeric.")
  }
  if (MinCorrelationScore < -1 | MinCorrelationScore > 1) {
    stop("MinCorrelationScore must be between -1 and 1.")
  }

  # Check the +/- window range
  if (!is.numeric(Window) | length(Window) != 1) {
    stop("Window must be a single numeric.")
  }
  Window <- abs(Window)

  ##################################
  ## Make the trelliscope display ##
  ##################################

  # Convert ProteoMatch class
  ProteoMatchTrelli <- Ms1Match
  class(ProteoMatchTrelli) <- c("data.table", "data.frame")

  # Filter ProteoMatch down to the correlation score
  ProteoMatchTrelli <- ProteoMatchTrelli %>%
    dplyr::filter(Correlation >= MinCorrelationScore)

  # List relevant ProteoMatch columns
  RelCol <- c("Protein", "Absolute Relative Error", "Correlation", "Molecular Formula",
              "Monoisotopic Mass", "Figure of Merit", "Charge", "Proteoform", "ID")

  # Calculate Median PPM Error and Minimum MZ
  MedianPPMError <- ProteoMatchTrelli %>%
    dplyr::select(ID, `PPM Error`, `M/Z`) %>%
    dplyr::group_by(ID) %>%
    dplyr::summarise(
      `Median PPM Error` = median(`PPM Error`),
      `Minimum Matched M/Z` = min(`M/Z`),
      `Peaks Matched` = length(`M/Z`)
    ) %>%
    dplyr::ungroup()
  
  # Generate trelliscope display
  ProteoMatchTrelli %>%
    dplyr::select(RelCol) %>%
    merge(MedianPPMError, by = "ID") %>%
    dplyr::mutate(ID = as.factor(ID)) %>%
    unique() %>%
    dplyr::mutate(
      panel = trelliscopejs::map_plot(ID, function(x) {plot_Ms1Match(PeakData, Ms1Match, x, Window, Trace = FALSE)})
    ) %>%
    trelliscopejs::trelliscope(
      path = Path,
      name = "MS1 Matches",
      nrow = 1,
      ncol = 1,
      thumb = T,
    )

}




