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
#' @param MatchingAlgorithm (character) Either "closest peak" or "highest abundance" where the "closest
#'    peak" implementation chooses the peak closest to the true M/Z value within the PPM window
#'    and "highest abundance" chooses the highest intensity peak within the PPM window. "closest peak"
#'    is recommended for peaks that have been peak picked with an external tool,
#'    and "highest abundance" is recommended for noisy datasets or those with many peaks.
#' @param MinAbundance (numeric) The minimum abundance (calculated intensity) permitted
#'     to be matched. Default is 0.1, which is 0.1%.
#' @param PPMThreshold (numeric) The window in PPM to search for peaks around the true peak. Required. Default is 10 ppm.
#' @param IsotopeRange (numeric) The minimum and maximum number of isotopes to consider. Default is c(3,20).
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
#'    PeakData = PeakData,
#'    MolecularFormula = MolForms_Test,
#'    MatchingAlgorithm = "closest peak",
#'    IsotopeRange = c(3,20)
#' )
#'
#' }
#'
#' @export
match_biomolecule_to_ms1 <- function(PeakData,
                                    MolecularFormulas,
                                    MatchingAlgorithm,
                                    MinAbundance = 0.1,
                                    PPMThreshold = 10,
                                    IsotopeRange = c(5, 20)) {

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

  # PPM Threshold
  if (!is.numeric(PPMThreshold) || PPMThreshold == 0) {
    stop("PPMThreshold must be a nonzero numeric.")
  }
  PPMThreshold <- abs(PPMThreshold)

  # MinAbundance should be a numeric value
  if (!is.numeric(MinAbundance) || MinAbundance < 0 | MinAbundance > 100) {
    stop("MinAbundance should be a numeric between 0 and 100, inclusive.")
  }

  # Check that max isotope is a numeric
  if (!is.numeric(IsotopeRange) | length(unique(abs(round(IsotopeRange)))) != 2) {
    stop("IsotopeRange must be a numeric with at least two unique values.")
  }
  IsotopeRange <- abs(round(IsotopeRange))
  MinIsotopes <- min(IsotopeRange)
  MaxIsotopes <- max(IsotopeRange)

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

    # Convert the molecular formula back to an atomic vector
    molform <- pspecterlib::as.molform(MolForm)

    # Get isotope profile (distribution). Match quality is determined by the limit function.
    IsoDist <- pspecterlib::calculate_iso_profile(molform, MinAbundance, limit = 0.001) %>%
      data.table::data.table() %>%
      dplyr::rename(`M/Z` = mass, Abundance = abundance, Isotope = isotope) %>%
      dplyr::select(-isolabel) %>%
      dplyr::mutate(`M/Z` = (`M/Z` + (Charge * AdductMass)) / Charge)

    # Add mass change if it is not NULL
    if (!is.null(MassShift) & MassShift != 0) {
      IsoDist$`M/Z` <- IsoDist$`M/Z` + (MassShift / Charge)
    }

    # Save the original IsoDist object before applying matches and filtering
    OrigIsoDist <- IsoDist

    ####################
    ## MATCH ISOTOPES ##
    ####################

    # Determine the theoretical mz tolerance
    IsoDist$`M/Z Search Window` <- IsoDist$`M/Z` * PPMThreshold / 1e6

    if (MatchingAlgorithm == "closest peak") {

      # For each theoretical peak, find the closest index in ms, where ms = theoretical
      LeftIndex <- findInterval(IsoDist$`M/Z`, PeakData$`M/Z`, rightmost.closed = FALSE, all.inside = TRUE)

      # Compute mz differences (absolute) to closest element to each side, smaller to the left and next greater to the right:
      IsoDist$`Left Difference` <- abs(PeakData$`M/Z`[LeftIndex] - IsoDist$`M/Z`)
      IsoDist$`Right Difference` <- abs(PeakData$`M/Z`[LeftIndex + 1] - IsoDist$`M/Z`)
      IsoDist$`Closest Index` <- LeftIndex

      # Set closest index as right side one, if difference is smaller:
      RightIndexBest <- which(IsoDist$`Right Difference` < IsoDist$`Left Difference`)
      IsoDist$`Closest Index`[RightIndexBest] <- IsoDist$`Closest Index`[RightIndexBest] + 1
      IsoDist$`M/Z Difference` <- abs(PeakData$`M/Z`[IsoDist$`Closest Index`] - IsoDist$`M/Z`)

      # Keep only matches within the tolerance
      IsoDist <- IsoDist[which(IsoDist$`M/Z Difference` < IsoDist$`M/Z Search Window`), ]
      IsoDist$`M/Z Experimental` <- PeakData$`M/Z`[IsoDist$`Closest Index`]
      IsoDist$`Intensity Experimental` <- PeakData$Intensity[IsoDist$`Closest Index`]
      IsoDist$`Abundance Experimental` <- PeakData$Abundance[IsoDist$`Closest Index`]

      # Remove non-necessary rows moving forward
      IsoDist <- IsoDist %>% dplyr::select(-c(`Left Difference`, `Right Difference`, `Closest Index`, `M/Z Difference`))

    } else if (MatchingAlgorithm == "highest abundance") {

      IsoDist <- IsoDist %>%
        dplyr::mutate(
          MZLower = `M/Z` - `M/Z Search Window`,
          MZUpper = `M/Z` + `M/Z Search Window`,
          `M/Z Experimental` = purrr::map2(MZLower, MZUpper, function(low, high) {
            sub <- PeakData[PeakData$`M/Z` >= low & PeakData$`M/Z` <= high,]
            if (nrow(sub) == 0) {return(NA)} else {return(sub[which.max(sub$Abundance), "M/Z"])}
          }) %>% unlist()
        ) %>%
        dplyr::filter(!is.na(`M/Z Experimental`)) %>%
        dplyr::mutate(
          `Intensity Experimental` = purrr::map(`M/Z Experimental`, function(x) {
            PeakData[PeakData$`M/Z` == x, "Intensity"]
          }) %>% unlist(),
          `Abundance Experimental` = purrr::map(`M/Z Experimental`, function(x) {
            PeakData[PeakData$`M/Z` == x, "Abundance"]
          }) %>% unlist()
        ) %>%
        dplyr::select(-c(MZLower, MZUpper))

    }

    # Calculate PPM Error
    IsoDist$`PPM Error` <- ((IsoDist$`M/Z Experimental` - IsoDist$`M/Z`) / IsoDist$`M/Z`) * 1e6

    # Subset down to numbers that are within 1 place of each other
    CloseValues <- IsoDist %>%
      dplyr::select(Isotope) %>%
      dplyr::mutate(
        Order = Isotope - dplyr::lag(Isotope),
        Order = ifelse(is.na(Order), 1, Order),
        Take = Order == dplyr::lag(Order) & Order == 1,
        Take = ifelse(is.na(Take), TRUE, Take),
        Take = ifelse(Take == TRUE, ifelse(dplyr::lag(Take) == FALSE, FALSE, TRUE), Take),
        Take = ifelse(is.na(Take), TRUE, Take)
      )

    IsoDist <- IsoDist[CloseValues$Take,]

    if (nrow(IsoDist) < min(IsotopeRange)) {return(NULL)}

    ##########################
    ## CALCULATED ABUNDANCE ##
    ##########################

    # Get max calculated intensity and max measured abundance
    calcInten <- unlist(IsoDist$Intensity)[which.max(unlist(IsoDist$Intensity))]
    maxAbun <- unlist(IsoDist$Abundance)[which.max(unlist(IsoDist$Abundance))]

    # Determine scale and scale intensity
    scalingFactor <- maxAbun / calcInten
    IsoDist$Abundance <- IsoDist$Intensity * scalingFactor

    ######################
    ## ABUNDANCE FILTER ##
    ######################

    # Flag abundance changes that are greater than the threshold
    IsoDist <- IsoDist %>%
      dplyr::mutate(
        `Abundance Diff` = `Abundance Experimental` -
          dplyr::lag(`Abundance Experimental`, default = dplyr::first(`Abundance Experimental`)),
        Flag = abs(`Abundance Diff`) >= MinAbundance | `Abundance Diff` == 0
      )

    # Determine where to subset from
    if (FALSE %in% IsoDist$Flag) {
      IsoDist <- IsoDist[1:(min(which(IsoDist$Flag == FALSE))-1),]
    }

    ######################
    ## CALCULATE SCORES ##
    ######################

    # Add correlation score
    if (nrow(IsoDist) >= min(IsotopeRange)) {

      # Generate an abundance match data.frame
      AbundanceDF <- merge(OrigIsoDist[,c("M/Z", "Abundance")], IsoDist[,c("M/Z", "Abundance Experimental")],
                           by = "M/Z", all.x = T)
      AbundanceDF$`Abundance Experimental`[is.na(AbundanceDF$`Abundance Experimental`)] <- 0

      # Calculate Absolute Relative Error, Correlation, and Figure of merit
      IsoDist$`Absolute Relative Error` <- 1/nrow(AbundanceDF) * sum(abs(AbundanceDF$Abundance - AbundanceDF$`Abundance Experimental`) / AbundanceDF$Abundance)
      IsoDist$Correlation <- stats::cor(AbundanceDF$`Abundance Experimental`, AbundanceDF$Abundance, method = "pearson")
      IsoDist$`Figure of Merit` <- nrow(AbundanceDF) / (sum((AbundanceDF$Abundance - AbundanceDF$`Abundance Experimental`)^2 + attributes(PeakData)$pspecter$MinimumAbundance^2))

      # Figure of merit
      IsoDist$`Figure of Merit` <- ifelse(is.infinite(IsoDist$`Figure of Merit`), NA, IsoDist$`Figure of Merit`)

      # Generate an identifier
      IsoDist$ID <- uuid::UUIDgenerate()

      # Add input columns
      IsoDist$`Monoisotopic Mass` <- MonoMass
      IsoDist$`Molecular Formula` <- MolForm
      IsoDist$Charge <- Charge
      IsoDist$`Mass Shift` <- MassShift
      IsoDist$AdductMasses <- AdductMass

      # Add missing columns and reorder
      return(IsoDist)

    } else {return(NULL)}

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
      AdductMasses = MolecularFormulas$AdductMasses[it]
    )
  }

  # If there is no MolFormTable, stop
  if (is.null(MolFormTable) || nrow(MolFormTable) == 0) {
    stop("No fragmentation patterns found. Try increasing the PPMThreshold, decreasing the minimum isotope range, lowering the noise filter, or lowering the minimum abundance.")
  }

  # Subset matches
  AllMatches <- merge(MolFormTable, MolecularFormulas, by = c("Monoisotopic Mass", "Mass Shift", "Molecular Formula", "Charge", "AdductMasses"), all.x = T) %>%
    dplyr::mutate(ID = as.numeric(as.factor(ID)))

  AllMatches <- AllMatches %>%
    dplyr::select(
      Identifiers, `M/Z`, `Monoisotopic Mass`, Abundance, Isotope, `M/Z Search Window`, `M/Z Experimental`,
      `Intensity Experimental`, `Abundance Experimental`, `PPM Error`, `Absolute Relative Error`,
      Correlation, `Figure of Merit`, Charge, Biomolecules, `Molecular Formula`,
      `Most Abundant Isotope`, ID
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