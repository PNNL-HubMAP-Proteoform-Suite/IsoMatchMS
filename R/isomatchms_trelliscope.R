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

#' Generate a trelliscope display for all Biomolecule to MS1 Matches
#'
#' @description [Trelliscope](https://github.com/hafen/trelliscopejs) allows
#'    users to filter and sort plots based on cognostic variables. For this dataset,
#'    the plots are generated with the plot_Ms1Match function, and the cognostics
#'    are: Identifiers, Absolute Relative Error, Correlation, Charge, Biomolecules, and ID.
#'
#' @param PeakData (peak_data object) A pspecterlib or data.table with "M/Z" and "Intensity". Required.
#' @param Ms1Match (IsoMatchMS_MatchedPeaks object) object from match_full_seq_ms1. Required.
#' @param Path (character) The base directory of the trelliscope application. Default is Downloads/Ms1Match.
#' @param MinCorrelationScore (numeric) The minimum correlation score to plot. Default is 0.7.
#' @param Window (numeric) The -/+ m/z value on either side of the matched spectra plot. Default is 2 m/z.
#'
#' @returns An html trelliscope display
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
#'   IsotopeMinimum = 2
#' )
#'
#' # Make the trelliscope display
#' isomatchms_trelliscope(PeakData = PeakData, Ms1Match = IsoMatch)
#'
#' }
#'
#'
#' @export
isomatchms_trelliscope <- function(PeakData,
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

  # Ms1Match should be a IsoMatchMS_MatchedPeaks object
  if ("IsoMatchMS_MatchedPeaks" %in% class(Ms1Match) == FALSE) {
    stop("Ms1Match must be a IsoMatchMS_MatchedPeaks object from match_proteoform_to_ms1")
  }

  # Check that minimum correlation score is within 0-1
  if (!is.na(MinCorrelationScore)) {
    if (!is.numeric(MinCorrelationScore)) {
      stop("MinCorrelationScore must be a numeric.")
    }
    if (MinCorrelationScore < -1 | MinCorrelationScore > 1) {
      stop("MinCorrelationScore must be between -1 and 1.")
    } 
  }

  # Check the +/- window range
  if (!is.numeric(Window) | length(Window) != 1) {
    stop("Window must be a single numeric.")
  }
  Window <- abs(Window)

  ##################################
  ## Make the trelliscope display ##
  ##################################

  # Convert IsoMatchMS class
  IsoMatchMSTrelli <- Ms1Match
  class(IsoMatchMSTrelli) <- c("data.table", "data.frame")
  IsoMatchMSTrelli <- IsoMatchMSTrelli %>% dplyr::rename(Identifier = Identifiers)

  # Filter IsoMatchMS down to the correlation score
  if (!is.na(MinCorrelationScore)) {
    IsoMatchMSTrelli <- IsoMatchMSTrelli %>%
      dplyr::filter(`Pearson Correlation` >= MinCorrelationScore) 
  }
  
  # List relevant IsoMatchMS columns
  RelCol <- c("Biomolecules", "Identifier", "Absolute Relative Error", "Pearson Correlation", "Molecular Formula",
              "Monoisotopic Mass", "Charge", "Adduct Name", "Adduct Mass", "Mass Shift", "ID")

  # Calculate Median PPM Error and Minimum MZ
  MedianPPMError <- IsoMatchMSTrelli %>%
    dplyr::select(ID, `PPM Error`, `M/Z`) %>%
    dplyr::group_by(ID) %>%
    dplyr::summarise(
      `Median PPM Error` = round(median(`PPM Error`, na.rm = T), 8),
      `Minimum Matched M/Z` = min(`M/Z`),
      `Peaks Matched` = sum(!is.na(`PPM Error`))
    ) %>%
    dplyr::ungroup()

  # Generate trelliscope display
  suppressWarnings({
    IsoMatchMSTrelli %>%
      dplyr::select(RelCol) %>%
      merge(MedianPPMError, by = "ID") %>%
      dplyr::mutate(ID = as.numeric(ID)) %>%
      unique() %>%
      dplyr::mutate(
        panel = trelliscopejs::map_plot(ID, function(x) {plot_Ms1Match(PeakData, Ms1Match, x, Window)})
      ) %>%
      dplyr::arrange(-`Pearson Correlation`) %>%
      trelliscopejs::trelliscope(
        path = Path,
        name = "MS1 Matches",
        nrow = 1,
        ncol = 1,
        thumb = T
      )
  })


}




