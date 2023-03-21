#' Matches experimental peaks to calculated isotope profile
#'
#' @description Returns the "IsoMatchMS_MatchedPeaks" object with peaks that
#'     have matched to the isotope profile. This object can be passed to downstream
#'     plotting functions.
#'
#' @param PeakData (peak_data object) from pspecterlib, made with "make_peak_data" or
#'     extracted from a raw or mzML file with "get_peak_data." Use of centroided data
#'     is recommended, but not required.
#' @param MolecularFormulas (IsoMatchMS_MolForm object) object with Molecular Formulas,
#'     Mass Shifts, and Charges.
#' @param AbundanceThreshold (numeric) The +/- percent abundance an isotope peak
#'     can vary and still be considered a match. For example, if the true isotope peak
#'     is 3 but the measured value is 4, an AbundanceThreshold of 50% will consider
#'     that peak a true match if it ranges between 1.5 - 4.5. Default is 50. 
#' @param PPMThreshold (numeric) The window in PPM to search for peaks around the true peak. Required. Default is 10 ppm.
#' @param IsotopeMinimum (numeric) The minimum number of isotopes to consider. Default is 3.
#' @param IsotopeAlgorithm (character) "isopat" uses the isopat package to calculate isotopes, 
#'     while "Rdisop" uses the Rdisop package. Though more accurate, Rdisop has been known 
#'     to crash on Windows computers when called iteratively more than 1000 times. 
#'     Default is Rdisop, though isopat is an alternative.
#'
#' @details
#' The data.table outputted by this function contains 12 columns
#' \tabular{ll}{
#' Identifier \tab The provided biomolecule identifier \cr
#' \tab \cr
#' Adduct \tab The provided adduct \cr
#' \tab \cr
#' M/Z \tab The calculated M/Z of a particular isotope \cr
#' \tab \cr
#' Intensity \tab The calculated relative intensity of a particular isotope \cr
#' \tab \cr
#' Isotope \tab The isotope number "n" in M+n \cr
#' \tab \cr
#' M/Z Search Window \tab Given the ppm, the M/Z Search Window is the window in which peaks are searched and matched \cr
#' \tab \cr
#' M/Z Experimental \tab The M/Z of the matched from the experimental peaks \cr
#' \tab \cr
#' Intensity Experimental \tab The intensity of the matched peak from the experimental peaks \cr
#' \tab \cr
#' PPM Error \tab How far off the calculated peak is from the experimental peak. Only peaks within the ppm threshold are kept \cr
#' \tab \cr
#' Absolute Relative Error \tab The sum of the absolute relative difference of calculated and experimental peak intensities. \cr
#' \tab \cr
#' Correlation \tab The cosine correlation of calculated and experimental peak intensities \cr
#' \tab \cr
#' Charge \tab The provided charges \cr
#' \tab \cr
#' Biomolecule \tab The provided biomolecule \cr
#' \tab \cr
#' ID \tab A unique ID for each Proteoform, Protein, and Charge combination used in plotting functions \cr
#' \tab \cr
#' }
#'
#' @importFrom foreach %dopar% foreach
#'
#' @returns (IsoMatchMS_MatchedPeaks object) A IsoMatchMS_MatchedPeaks object, which is a data.table containing the
#'     Identifier, Adduct, M/Z, Intensity, Isotope, M/Z Search Window, Experimental M/Z,
#'     Experimental Intensity, PPM Error, Absolute Relative Error, Correlation,
#'     Charge, Biomolecule, and an ID.
#'
#' @examples
#' \dontrun{
#'
#' # Run two examples with two charge states
#' MolForms_Test <- calculate_molform(
#'    Biomolecules = c("M.SS[Methyl]S.V", "M.SS[6]S[7].V"),
#'    BioType = "ProForma",
#'    Identifiers = c("Test1", "Test2"),
#'    Charge = 1:2
#' )
#'
#' # Generate some experimental peak data to match
#' PeakData <- pspecterlib::make_peak_data(
#'    MZ = c(294.1296, 295.1325, 296.1343, 297.1369, 298.1390),
#'    Intensity = c(868.3680036, 110.9431876, 18.7179196, 1.7871629, 0.1701294)
#' )
#'
#' # Run algorithm
#' match_biomolecule_to_ms1(
#'     PeakData = PeakData,
#'     MolecularFormula = MolForms_Test,
#'     IsotopeMinimum = 2
#' )
#'
#' }
#'
#' @export
match_biomolecule_to_ms1 <- function(PeakData,
                                     MolecularFormulas,
                                     AbundanceThreshold = 50,
                                     PPMThreshold = 10,
                                     IsotopeMinimum = 3,
                                     IsotopeAlgorithm = "Rdisop") {

  ##################
  ## CHECK INPUTS ##
  ##################

  # Check that peak data is of the appropriate class
  if ("peak_data" %in% class(PeakData) == FALSE) {
    stop("PeakData must be a pspecterlib peak_data object.")
  }

  # Check that molecular formula is a string
  if (inherits(MolecularFormulas, "IsoMatchMS_MolForm") == FALSE) {
    stop("MolecularFormula must be an IsoMatchMS_MolForm object.")
  }

  # Abundance Threshold 
  if (!is.numeric(AbundanceThreshold) || AbundanceThreshold == 0) {
    stop("AbundanceThreshold must be a nonzero numeric.")
  }
  AbundanceThreshold <- abs(AbundanceThreshold)
  
  # PPM Threshold
  if (!is.numeric(PPMThreshold) || PPMThreshold == 0) {
    stop("PPMThreshold must be a nonzero numeric.")
  }
  PPMThreshold <- abs(PPMThreshold)

  # Check that minimum isotope is a numeric
  if (!is.numeric(IsotopeMinimum) && IsotopeMinimum != 0) {
    stop("IsotopeRange must be a non-zero numeric value.")
  }
  IsotopeRange <- abs(round(IsotopeMinimum))
  
  ##################
  ## RUN ITERATOR ##
  ##################

  .match_proteoform_to_ms1_iterator <- function(MonoMass,
                                                MolForm,
                                                Charge,
                                                MassShift,
                                                AdductMass) {

    ########################
    ## CALCULATE ISOTOPES ##
    ########################

    # Make sure that the PeakData overlaps with +10 the monoisotopic mass
    if (nrow(PeakData[PeakData$`M/Z` >= round(MonoMass) & PeakData$`M/Z` <= round(MonoMass)+10,]) == 0) {
      return(NULL)
    }
    
    # Get isotope profile (distribution). Match quality is determined by the limit function.
    IsoDist <- pspecterlib::calculate_iso_profile(
      molform = pspecterlib::as.molform(MolForm),
      algorithm = IsotopeAlgorithm
    ) %>%
      dplyr::rename(`M/Z` = mass, Abundance = abundance, Isotope = isotope) %>%
      dplyr::select(-isolabel) %>%
      dplyr::mutate(`M/Z` = (`M/Z` + (Charge * AdductMass)) / Charge)
    
    # If we don't calculate the minimum number of peaks, dispose 
    if (nrow(IsoDist) < min(IsotopeRange)) {return(NULL)}

    # Add mass change if it is not NULL
    if (!is.null(MassShift) & MassShift != 0) {
      IsoDist$`M/Z` <- IsoDist$`M/Z` + (MassShift / Charge)
    }

    # Save the original IsoDist object before applying matches and filtering
    OrigIsoDist <- IsoDist

    ####################
    ## MATCH ISOTOPES ##
    ####################

    # Determine the theoretical mz tolerance for each value
    IsoDist$`M/Z Search Window` <- IsoDist$`M/Z` * PPMThreshold / 1e6

    # Set the M/Z search window
    IsoDist <- IsoDist %>%
      dplyr::mutate(
        MZLower = `M/Z` - `M/Z Search Window`,
        MZUpper = `M/Z` + `M/Z Search Window`
      )
    
    # Subset PeakData
    PeakSub <- PeakData[PeakData$`M/Z` >= min(IsoDist$MZLower) & PeakData$`M/Z` <= max(IsoDist$MZUpper),]
    
    # Rescale abundances on the PeakData
    PeakRe <- pspecterlib::make_peak_data(MZ = PeakSub$`M/Z`, Intensity = PeakSub$Abundance) 
    class(PeakRe) <- c("data.table", "data.frame")
    
    # Match the experimental peaks based on closest MZ and Abundance 
    IsoDist <- IsoDist %>% dplyr::mutate(
      
      `Closest Index` = purrr::pmap(list(MZLower, MZUpper, Abundance), function(low, high, abun) {
        
        # Get peak range 
        sub <- PeakRe[PeakRe$`M/Z` >= low & PeakRe$`M/Z` <= high & PeakRe$Abundance >= abun - (abun * AbundanceThreshold/100) & PeakRe$Abun <= abun + (abun * AbundanceThreshold/100) ,]
        
        # If no peaks, return 0
        if (nrow(sub) == 0) {return(NA)}
        
        # Make rank table and identify best match. It should be the closest MZ and Abundance,
        # with a preference for MZ 
        BestMatch <- data.frame(
          Names = 1:nrow(sub),
          Rank1 = order(abs(sub$`M/Z` - mean(c(low, high)))),
          Rank2 = order(abs(sub$Abundance - abun))
        ) %>% 
          dplyr::mutate(
            TotalRank = Rank1 + Rank2
          ) %>%
          dplyr::filter(TotalRank == min(TotalRank)) %>%
          dplyr::filter(Rank1 == min(Rank1)) %>% 
          dplyr::select(Names) %>%
          unlist()
        
        return(which.min(abs(PeakRe$`M/Z` - sub$`M/Z`[BestMatch])))
        
      }) %>% unlist(),
      
      `M/Z Experimental` = ifelse(!is.na(`Closest Index`), PeakSub$`M/Z`[`Closest Index`], NA),
      `Intensity Experimental` = ifelse(!is.na(`Closest Index`), PeakSub$Intensity[`Closest Index`], NA),
      `Abundance Experimental` =  ifelse(!is.na(`Closest Index`), PeakSub$Abundance[`Closest Index`], NA),
      `PPM Error` = ifelse(!is.na(`Closest Index`), (`M/Z Experimental` - `M/Z`) / `M/Z` * 1e6, NA)
  
    ) %>%
      dplyr::select(-`Closest Index`)
    
    # Make sure we have enough matches 
    if ((IsoDist %>% dplyr::filter(!is.na(`M/Z Experimental`)) %>% nrow()) < IsotopeMinimum) {
      return(NULL)
    }
        
    ######################
    ## CALCULATE SCORES ##
    ######################

    # Generate an abundance match data.frame
    AbundanceDF <- merge(OrigIsoDist[,c("M/Z", "Abundance")], IsoDist[,c("M/Z", "Abundance Experimental")], by = "M/Z", all.x = T)
    AbundanceDF$`Abundance Experimental`[is.na(AbundanceDF$`Abundance Experimental`)] <- 0

    # Calculate Absolute Relative Error and Pearson Correlation
    IsoDist$`Absolute Relative Error` <- round(1/nrow(AbundanceDF) * sum(abs(AbundanceDF$Abundance - AbundanceDF$`Abundance Experimental`) / AbundanceDF$Abundance), 8)
    IsoDist$`Pearson Correlation` <- round(stats::cor(AbundanceDF$`Abundance Experimental`, AbundanceDF$Abundance, method = "pearson"), 8)

    # Generate an identifier
    IsoDist$ID <- uuid::UUIDgenerate()

    # Add input columns
    IsoDist$`Monoisotopic Mass` <- MonoMass
    IsoDist$`Molecular Formula` <- MolForm
    IsoDist$Charge <- Charge
    IsoDist$`Mass Shift` <- MassShift
    IsoDist$`Adduct Mass` <- AdductMass

    # Add missing columns and reorder
    return(IsoDist)

  }

  # Implement parallel computing for speed
  doParallel::registerDoParallel(parallel::detectCores())
  
  # Iterate through and match molecular formula data. Remove NULLs and pull out the ID
  MolFormTable <- foreach(it = 1:nrow(MolecularFormulas), .combine = rbind) %dopar% {
    .match_proteoform_to_ms1_iterator(
      MonoMass = MolecularFormulas$`Monoisotopic Mass`[it],
      MolForm = MolecularFormulas$`Molecular Formula`[it],
      Charge = MolecularFormulas$Charge[it],
      MassShift = MolecularFormulas$`Mass Shift`[it],
      AdductMass = MolecularFormulas$`Adduct Mass`[it]
    )
  }
  
  # Saving this chunk for debugging
  
  #MolFormTable <- do.call(rbind, lapply(1:nrow(MolecularFormulas), function(it) {
  #  message(it)
  #  .match_proteoform_to_ms1_iterator(
  #      MonoMass = MolecularFormulas$`Monoisotopic Mass`[it],
  #      MolForm = MolecularFormulas$`Molecular Formula`[it],
  #      Charge = MolecularFormulas$Charge[it],
  #      MassShift = MolecularFormulas$`Mass Shift`[it],
  #      AdductMass = MolecularFormulas$Adduct[it]
  #  )
  #}))

  # If there is no MolFormTable, stop
  if (is.null(MolFormTable) || nrow(MolFormTable) == 0) {
    stop("No fragmentation patterns found. Try increasing the PPMThreshold, decreasing the minimum isotope range, lowering the noise filter, or lowering the minimum abundance.")
  }

  # Subset matches
  AllMatches <- merge(MolFormTable, 
                      MolecularFormulas, 
                      by = c("Monoisotopic Mass", "Mass Shift", "Molecular Formula", "Charge", "Adduct Mass"), 
                      all.x = T) %>%
    dplyr::mutate(ID = as.numeric(as.factor(ID)))

  AllMatches <- AllMatches %>%
    dplyr::select(
      Identifiers, `Adduct Mass`, `Adduct Name`, `M/Z`, `Mass Shift`, `Monoisotopic Mass`, 
      Abundance, Isotope, `M/Z Search Window`, `M/Z Experimental`,
      `Intensity Experimental`, `Abundance Experimental`, `PPM Error`, 
      `Absolute Relative Error`, `Pearson Correlation`, Charge, Biomolecules, 
      `Molecular Formula`, `Most Abundant Isotope`, ID
     ) %>%
    unique()

  # If no matches, return NULL
  if (is.null(AllMatches)) {
    message("No matches detected.")
    return(NULL)
  }

  # Add class, and return object
  class(AllMatches) <- c(class(AllMatches), "IsoMatchMS_MatchedPeaks")

  return(AllMatches)

}
