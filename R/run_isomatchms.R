#' Run the IsoMatchMS pipeline with a list of PTMs, summed peak data, and algorithm settings
#'
#' @description A wrapper for quickly running all the steps of the IsoMatchMS pipeline
#'
#' @param Biomolecules A vector of Molecular Formulas or ProForma strings (i.e. "M.AA[Acetyl⁠]⁠AA[3.2⁠]⁠.V").
#'    ProForma strings can be pulled from mzid files (MS-GF+, MSPathFinder) with pull_modifications_from_mzid, 
#'    or made with create_proforma for MSPathFinder, ProSight, and pTop modifications. TopPIC proteoforms 
#'    are provided as ProForma strings. Required.
#' @param BioType A string indicating whether the Biomolecules are "ProForma" strings or "Molecular Formula". Required. 
#' @param SummedSpectra A pspecterlib peak_data object of all summed MS1 spectra. This
#'    can be generated with an mzML file and the sum_ms1_spectra function, pulled from
#'    a summed mzML with pspecterlib::get_peak_data, or simply created with 
#'    pspecterlib::make_peak_data. Required.  
#' @param SettingsFile Path to a xlsx file with all user-set parameters. For examples, 
#'    try xlsx::read.xlsx(system.file("extdata", "Intact_Protein_Defaults.xlsx", package = "IsoMatchMS"), 1). 
#'    Required.
#' @param Path The base directory of the output trelliscope application. Default is Downloads/Ms1Match.
#' @param Identifiers A vector of identifiers for each biomolecule (i.e. protein, glycan, etc.). Optional.
#' @param Messages A TRUE/FALSE indicating whether messages should be printed. 
#'
#' @details Before running the main algorithm, you will need Molecular Formulas (or ProForma strings)
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
#' Documentation for example intact proteomic data can be found [here](https://www.sciencedirect.com/science/article/pii/S1535947622002997)
#' 
#' Documentation for example digested proteomic data can be found [here](https://google.com)
#' 
#' Documentation for example glycan data can be found [here](https://google.com)
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
#' ####################
#' ## SHORT EXAMPLES ##
#' ####################
#' 
#' ## INTACT PROTEOMICS EXAMPLE ##
#' 
#' protein_data <- read.csv(system.file("extdata", "Intact_Proteins_List_Short.csv", package = "IsoMatchMS"))
#' peak_data <- readRDS(system.file("extdata", "Intact_PeakData.RDS", package = "IsoMatchMS"))
#' 
#' 
#' ## DIGESTED PROTEOMICS EXAMPLE ##
#' 
#' ###################
#' ## FULL EXAMPLES ##
#' ###################
#' 
#' ## INTACT PROTEOMICS EXAMPLE ##
#' 
#' # Extract the summed peak data with any tool (here, pspecterlib, which uses mzR)
#' MSPath <- system.file("extdata", "Intact_Protein_Summed_MS1.mzML", package = "IsoMatchMS")
#' ScanMetadata <- pspecterlib::get_scan_metadata(MSPath)
#' PeakData <- pspecterlib::get_peak_data(ScanMetadata, ScanNumber = 1, MinAbundance = 0.1)
#' 
#' # Now, get the molecular formulas (required) and protein names (optional)
#' ProteinData <- read.csv(system.file("extdata", "Intact_Proteins_List.csv", package = "IsoMatchMS"), header = T)
#' Proteoform <- ProteinData$Proteoform
#' ProteinNames <- ProteinData$Protein.accession
#' 
#' # Get the settings file path
#' Intact_Settings_Path <- system.file("extdata", "Intact_Protein_Defaults.xlsx", package = "IsoMatchMS")
#' 
#' # Now, run the main pipeline 
#' run_isomatchms(
#'     Biomolecule = Proteoform,
#'     BioType = "ProForma",
#'     SummedSpectra = PeakData, 
#'     SettingsFile = Intact_Settings_Path, 
#'     Path = file.path(.getDownloadsFolder(), "Intact_Example"),
#'     Identifiers = ProteinNames,
#'     Messages = TRUE
#' )
#' 
#' ## DIGESTED PROTEOMICS EXAMPLE ##
#' 
#' # Pull summed peak data from CSV
#' Peaks <- read.csv(system.file("extdata", "Peptides_PeakData.csv", package = "IsoMatchMS"))
#' PeakData <- pspecterlib::make_peak_data(MZ = Peaks$M.Z, Intensity = Peaks$Intensity)
#' 
#' # Now, get the molecular formulas (required) and protein names (optional)
#' ProteinData <- read.csv(system.file("extdata", "Peptides_List.csv", package = "IsoMatchMS"))
#' ProForma <- ProteinData$Proteoform
#' ProteinNames <- ProteinData$Protein
#' 
#' # Get the settings file path
#' Digested_Settings_Path <- system.file("extdata", "Peptides_Defaults.xlsx", package = "IsoMatchMS")
#' 
#' # Now, run the main pipeline 
#' run_isomatchms(
#'     Biomolecule = ProForma,
#'     BioType = "ProForma",
#'     SummedSpectra = PeakData, 
#'     SettingsFile = Digested_Settings_Path, 
#'     Path = file.path(.getDownloadsFolder(), "Digested_Example"),
#'     Identifiers = ProteinNames,
#'     Messages = TRUE
#' )
#' 
#' }
#'
#' @export
run_isomatchms <- function(Biomolecules,
                            BioType,
                            SummedSpectra,
                            SettingsFile,
                            Path = file.path(.getDownloadsFolder(), "Ms1Match"),
                            Identifiers = NULL,
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
    Biomolecules = Biomolecules,
    BioType = BioType, 
    Identifiers = Identifiers,
    Charge = Settings[Settings$Parameter == "Charges", "Default"] %>% strsplit(",") %>% unlist() %>% as.numeric(),
    AddMostAbundantIsotope = Settings[Settings$Parameter == "AddMAI", "Default"] %>% unlist() %>% as.logical(),
    AdductMasses = Settings[Settings$Parameter == "ProtonMass", "Default"] %>% as.numeric()
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
  MatchedPeaks <- match_biomolecule_to_ms1(
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
  isomatchms_trelliscope(
    PeakData = FilteredData,
    Ms1Match = MatchedPeaks,
    Path = file.path(Path, "Trelliscope"),
    MinCorrelationScore = Settings[Settings$Parameter == "CorrelationMinimum", "Default"] %>% as.numeric(),
    Window = Settings[Settings$Parameter == "PlottingWindow", "Default"] %>% as.numeric()
  )

}
