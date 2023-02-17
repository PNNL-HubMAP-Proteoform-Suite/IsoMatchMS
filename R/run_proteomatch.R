#' Run the ProteoMatch pipeline with a list of PTMs, summed peak data, and algorithm settings
#'
#' @description A wrapper for quickly running all the steps of the ProteoMatch pipeline
#'
#' @param Modifications A vector of ProForma strings (i.e. "M.AA`[`Acetyl`]`AA`[`3.2`]`.V") 
#'    or Molecular Formulas. ProForma strings can be pulled from mzid files (MS-GF+, MSPathFinder) 
#'    with pull_modifications_from_mzid, or made with create_proforma for MSPathFinder,
#'    ProSight, and pTop modifications. TopPIC proteoforms are provided as ProForma
#'    strings. Required.
#' @param ModType A string indicating whether the Modifications are "ProForma" strings or "Molecular Formula". Required. 
#' @param SummedSpectra A pspecterlib peak_data object of all summed MS1 spectra. This
#'    can be generated with an mzML file and the sum_ms1_spectra function, pulled from
#'    a summed mzML with pspecterlib::get_peak_data, or simply created with 
#'    pspecterlib::make_peak_data. Required.  
#' @param SettingsFile Path to a xlsx file with all user-set parameters. For examples, 
#'    try xlsx::read.xlsx(system.file("extdata", "Intact_Protein_Defaults.xlsx", package = "ProteoMatch"), 1). 
#'    Required.
#' @param Path The base directory of the output trelliscope application. Default is Downloads/Ms1Match.
#' @param Proteins An optional list of protein names. Used in the trelliscope display. Must
#'    be the same length and in the same order as Modifications. Default is NULL.
#' @param Messages A TRUE/FALSE indicating whether messages should be printed. 
#'
#' @details Before running the main algorithm, you will need ProForma strings (or molecular formulas)
#' and a summed spectra. There are several ways to get ProForma strings, including
#' pull_modifications_from_mzid and create_proforma. Spectra may be summed with sum_ms1_spectra,
#' or uploaded from an mzML (see pspecterlib::get_scan_metadata and pspecterlib::get_peak_data) or
#' created from a csv (see pspecterlib::make_peak_data). 
#' 
#' The algorithm is run in 4 main steps:
#' 1. Molecular Formulas and masses are calculated with calculate_molform
#' 2. Input summed spectra are filtered with filter_peaks
#' 3. Isotopic distributions are calculated with the molecular formulas, and matched
#'     to the filtered spectra
#' 4. All high scoring matches are visualized with a trelliscope display
#' 
#' Intact proteomic data can be found [here](https://www.sciencedirect.com/science/article/pii/S1535947622002997)
#' 
#' Digested proteomic data can be found [here](https://google.com)
#' 
#' @returns 
#' \tabular{ll}{
#' MolecularFormulas \tab A CSV with all calculated molecular formulas from the ProForma strings. \cr
#' \tab \cr
#' FilteredPeaks \tab A CSV with the final spectrum used after filtering \cr 
#' \tab \cr
#' Matched_Isotope_Distributions \tab A CSV with all matched isotope distributions, regardless of quality \cr 
#' \tab \cr
#' Trelliscope_Display \tab Visualization of all matched isotope distributions above the quality threshold \cr
#' \tab \cr
#' }
#' 
#' @examples 
#' \dontrun{
#' 
#' ## INTACT PROTEOMICS EXAMPLE ##
#' 
#' # Extract the summed peak data with any tool (here, pspecterlib, which uses mzR)
#' MSPath <- system.file("extdata", "Intact_Protein_Summed_MS1.mzML", package = "ProteoMatch")
#' ScanMetadata <- pspecterlib::get_scan_metadata(MSPath)
#' PeakData <- pspecterlib::get_peak_data(ScanMetadata, ScanNumber = 1, MinAbundance = 0.1)
#' 
#' # Now, get the molecular formulas (required) and protein names (optional)
#' ProteinData <- read.table(system.file("extdata", "Intact_Protein_Molecular_Formulas.tsv", package = "ProteoMatch"), sep = "\t", header = T)
#' MolecularFormulas <- ProteinData$Formula
#' ProteinNames <- ProteinData$Protein
#' 
#' # Get the settings file path
#' Intact_Settings_Path <- system.file("extdata", "Intact_Protein_Defaults.xlsx", package = "ProteoMatch")
#' 
#' # Now, run the main pipeline 
#' run_proteomatch(
#'     Modifications = MolecularFormulas,
#'     ModType = "Molecular Formula",
#'     SummedSpectra = PeakData, 
#'     SettingsFile = Intact_Settings_Path, 
#'     Path = file.path(.getDownloadsFolder(), "Intact_ProteoMatch_Example"),
#'     Proteins = ProteinNames,
#'     Messages = TRUE
#' )
#' 
#' ## DIGESTED PROTEOMICS EXAMPLE ##
#' 
#' # Pull summed peak data from CSV
#' Peaks <- read.csv(system.file("extdata", "Digested_Protein_PeakData.csv", package = "ProteoMatch"))
#' PeakData <- pspecterlib::make_peak_data(MZ = Peaks$M.Z, Intensity = Peaks$Intensity)
#' 
#' # Now, get the molecular formulas (required) and protein names (optional)
#' ProteinData <- read.csv(system.file("extdata", "Digested_Proteins.csv", package = "ProteoMatch"))
#' ProForma <- ProteinData$Proteoform
#' ProteinNames <- ProteinData$Protein
#' 
#' # Get the settings file path
#' Digested_Settings_Path <- system.file("extdata", "Digested_Defaults.xlsx", package = "ProteoMatch")
#' 
#' # Now, run the main pipeline 
#' run_proteomatch(
#'     Modifications = ProForma,
#'     ModType = "ProForma",
#'     SummedSpectra = PeakData, 
#'     SettingsFile = Digested_Settings_Path, 
#'     Path = file.path(.getDownloadsFolder(), "Digested_ProteoMatch_Example"),
#'     Proteins = ProteinNames,
#'     Messages = TRUE
#' )
#' 
#' }
#'
#' @export
run_proteomatch <- function(Modifications,
                            ModType,
                            SummedSpectra,
                            SettingsFile,
                            Path = file.path(.getDownloadsFolder(), "Ms1Match"),
                            Proteins = NULL,
                            Messages = TRUE) {

  ##################
  ## CHECK INPUTS ##
  ##################
  
  # Modifications, ModType, and Proteins are checked by the calculate_molform function
  
  # SummedSpectra are checked in the filter_peaks function

  # The settings file should have be an xlsx
  if (!grepl(".xlsx", SettingsFile)) {
    stop("SettingsFile needs to be an xlsx file.")
  }
  Settings <- suppressWarnings({xlsx::read.xlsx(SettingsFile, 1)})

  # The settings file should have all the required parameters
  RequiredRow <- c("MZRange", "NoiseFilter", "Charges", "MatchingAlgorithm", "MinimumAbundance", "CorrelationMinimum",
    "PPMThreshold", "AddMAI", "IsotopeRange", "PlottingWindow", "ProtonMass")
  if (!all(Settings$Parameter %in% RequiredRow)) {
    stop("Settings file is missing: ",
         paste0(RequiredRow[!RequiredRow  %in% Settings$Parameter], ", ", collapse = ""))
  }
  
  # Messages should be a TRUE or FALSE
  if (!is.logical(Messages) && (Messages != TRUE | Messages != FALSE)) {
    stop("Messages should be a TRUE or FALSE.")
  }

  ##################
  ## RUN PIPELINE ##
  ##################

  # 0. Create output directory
  if (!dir.exists(Path)) {dir.create(Path)}
  
  # 1. Calculate Molecular Formula
  if (Messages) {message("Calculating molecular formulas...")}
  MolForm <- calculate_molform(
    Modifications = Modifications,
    ModType = ModType, 
    Proteins = Proteins,
    Charge = Settings[Settings$Parameter == "Charges", "Default"] %>% strsplit(",") %>% unlist() %>% as.numeric(),
    AddMostAbundantIsotope = Settings[Settings$Parameter == "AddMAI", "Default"] %>% unlist() %>% as.logical(),
    ProtonMass = Settings[Settings$Parameter == "ProtonMass", "Default"] %>% as.numeric()
  )
  write.csv(MolForm, file.path(Path, "Molecular_Formulas.csv"), row.names = F, quote = F)

  # 2. Filter the data
  if (Messages) {message("Filtering peaks...")}
  FilteredData <- filter_peaks(
    PeakData = SummedSpectra,
    MZRange = Settings[Settings$Parameter == "MZRange", "Default"] %>% strsplit("-") %>% unlist() %>% as.numeric(),
    NoiseFilter = Settings[Settings$Parameter == "NoiseFilter", "Default"] %>% as.numeric()
  )
  write.csv(FilteredData, file.path(Path, "Filtered_Peaks.csv"), row.names = F, quote = F)
  
  # 3. Match Peaks
  if (Messages) {message("Matching spectra...")}
  MatchedPeaks <- match_proteoform_to_ms1(
    PeakData = FilteredData,
    MolecularFormulas = MolForm,
    MatchingAlgorithm = Settings[Settings$Parameter == "MatchingAlgorithm", "Default"] %>% unlist(),
    MinAbundance = Settings[Settings$Parameter == "MinimumAbundance", "Default"] %>% unlist() %>% as.numeric(),
    PPMThreshold = Settings[Settings$Parameter == "PPMThreshold", "Default"] %>% as.numeric(),
    IsotopeRange = Settings[Settings$Parameter == "IsotopeRange", "Default"] %>% strsplit(",") %>% unlist() %>% as.numeric(),
    ProtonMass = Settings[Settings$Parameter == "ProtonMass", "Default"] %>% as.numeric()
  )
  if (is.null(MatchedPeaks)) {
    write.csv("No matches found", file.path(Path, "Matched_Isotope_Distributions.csv"), row.names = F, quote = F)
    return(NULL)
  }
  write.csv(MatchedPeaks, file.path(Path, "Matched_Isotope_Distributions.csv"), row.names = F, quote = F)

  # 4. Make the trelliscope display
  if (Messages) {message("Generating trelliscope...")}
  proteomatch_trelliscope(
    PeakData = FilteredData,
    Ms1Match = MatchedPeaks,
    Path = file.path(Path, "Trelliscope"),
    MinCorrelationScore = Settings[Settings$Parameter == "CorrelationMinimum", "Default"] %>% as.numeric(),
    Window = Settings[Settings$Parameter == "PlottingWindow", "Default"] %>% as.numeric()
  )

}
