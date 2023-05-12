#' Generate a Molecular Formula Table from provided molecular formulas, or ProForma strings
#'
#' @description Returns the "IsoMatchMS_MolForm" object with each biomolecule's molecular formula and adjusted mass.
#'
#' @param Biomolecules (character) A vector of Molecular Formulas or ProForma strings (i.e. "M.AA`[`Acetyl`]`AA`[`3.2`]`.V").
#'    ProForma strings can be pulled from mzid files (MS-GF+, MSPathFinder)
#'    with pull_modifications_from_mzid, or made with create_proforma for MSPathFinder,
#'    ProSight, and pTop modifications. TopPIC proteoforms are provided as ProForma
#'    strings. Required.
#' @param BioType (character) A string indicating whether the Biomolecules are "ProForma" strings or "Molecular Formula". Required.
#' @param Identifiers (character) A vector of identifiers for each biomolecule (i.e. protein, glycan, etc.). Optional.
#' @param Charge (numeric) The range of charges to test. Default is 1.
#' @param AddMostAbundantIsotope (boolean) A flag to determine whether the Most Abundant Isotope (MAI) should be calculated for
#'     every function. This parameter will slow down tool. Default is FALSE.
#' @param AdductMasses (vector) A named vector of the masses of adducts to be tested. 
#'     A max of 5 masses can be given. Proton Adducts are the default. Default is c(proton = 1.00727647).
#' @param IsotopeAlgorithm (character) "isopat" uses the isopat package to calculate 
#'     isotopes, while "Rdisop" uses the Rdisop package. Though more accurate, 
#'     Rdisop has been known to crash on Windows computers when called iteratively 
#'     more than 1000 times. Default is Rdisop, though isopat is an alternative.
#' @param MinAbundance (numeric) The minimum abundance (calculated intensity) permitted to be matched.
#'     Default is 0.1, which is 0.1%. Used for most abundant isotope. This is a pspecterlib-specific
#'     parameter and shouldn't need to be changed for IsoMatchMS.
#'
#' @details
#' The data.table outputted by this function returns 8 columns.
#' \tabular{ll}{
#' Biomolecules \tab The provided biomolecule molecular formulas or proforma strings \cr
#' \tab \cr
#' Identifiers \tab The provided biomolecule identifier \cr
#' \tab \cr
#' Adduct Names \tab The provided adduct names \cr
#' \tab \cr
#' Adduct Masses \tab The provided adduct masses \cr
#' \tab \cr
#' Charge \tab The provided charges \cr
#' \tab \cr
#' Molecular Formula \tab The molecular formula of the biomolecule, or derived from the peptide/protein sequence and its modifications \cr
#' \tab \cr
#' Mass Shift \tab The total mass change indicated in the proforma string, if provided \cr
#' \tab \cr
#' Monoisotopic Mass \tab The biomolecule's monoisotopic peak based on molecular formula and mass shift \cr
#' \tab \cr
#' Most Abundant Isotope \tab The biomolecule's most abundant isotope based on molecular formula and mass shift. Only added if requested. It will be calculated during the isotope matching. \cr
#' \tab \cr
#' }
#'
#' @importFrom foreach %dopar% foreach
#'
#' @returns (IsoMatchMS_MolForm) A IsoMatchMS_MolForm object, which is a data.table containing the
#'     molecular formula, mass shift, monoisotopic mass, most abundant isotope,
#'     biomolecule identifier, and charge.
#'
#' @examples
#' \dontrun{
#'
#' # Run one example with three charge states
#' calculate_molform(
#'    Biomolecules = "M.(S)[Acetyl]ATNNIAQARKLVEQLRIEAGIERIKVSKAASDLMSYCEQHARNDPLLVGVPASENPFKDK(KPCIIL)[-52.9879].",
#'    BioType = "ProForma",
#'    Identifiers = "O60262",
#'    Charge = 1:3,
#'    AddMostAbundantIsotope = TRUE
#' )
#'
#' # Run two examples with two charge states
#' calculate_molform(
#'    Biomolecules = c("M.SS[Methyl]S.V", "M.S[Methyl]S[22]S[23].V"),
#'    BioType = "ProForma",
#'    Charge = 1:2,
#' )
#'
#' # Run an example with molecular formulas and an adduct
#' calculate_molform(
#'    Biomolecules = c("C6H12O6", "C2H4O1"),
#'    BioType = "Molecular Formula",
#'    Identifiers = c("Glucose", "Acetyl"),
#'    AdductMasses = c(proton = 1.00727637, sodium = 22.989769),
#'    AddMostAbundantIsotope = TRUE
#' )
#'
#' }
#'
#' @export
calculate_molform <- function(Biomolecules,
                              BioType,
                              Identifiers = NULL,
                              Charge = 1,
                              AddMostAbundantIsotope = FALSE,
                              AdductMasses = c(proton = 1.00727647),
                              IsotopeAlgorithm = "Rdisop",
                              MinAbundance = 0.1) {

  ##################
  ## CHECK INPUTS ##
  ##################

  # Modifications should be a string
  if (is.null(Biomolecules) || !is.character(Biomolecules)) {
    stop("Biomolecules must be a vector of characters.")
  }

  # Blank or empty modifications are not permitted
  if (lapply(Biomolecules, function(x) {is.null(x) || is.na(x) || x == ""}) %>% unlist() %>% any()) {
    stop("Biomolecules cannot be NULL, NA, or blank.")
  }

  # Check that BioType is "Molecular Formula" or "ProForma"
  if (length(BioType) != 1 | !is.character(BioType) || (BioType != "ProForma" & BioType != "Molecular Formula")) {
    stop("BioType can either be 'ProForma' or 'Molecular Formula'.")
  }

  # If Identifiers is not NULL...
  if (!is.null(Identifiers)) {

    # Identifiers should be a string
    if (!is.character(Identifiers)) {
      stop("Identifiers must be a vector of characters.")
    }

    # Biomolecules and Identifiers should be the same length
    if (length(Biomolecules) != length(Identifiers)) {
      stop("Biomolecules and Identifiers must be of the same length.")
    }

  } else {
    Identifiers <- rep("NA", length(Biomolecules))
  }

  # Charge must be numeric and will be rounded to integers
  if (!is.numeric(Charge)) {
    stop("Charge should be numeric.")
  }
  Charge <- unique(round(abs(Charge)))

  # Check that AddMostAbundantIsotope is a TRUE/FALSE
  if (!is.logical(AddMostAbundantIsotope) | is.na(AddMostAbundantIsotope)) {
    stop("AddMostAbundnantIsotope must be true or false.")
  }

  # AdductMasses must be a single numeric
  if (length(AdductMasses) < 1 | length(AdductMasses) > 5) {
    stop("AdductMasses must be a list of length 1 to 5.")
  }

  # AdductMasses must be must contain numerics
  for (x in AdductMasses) {
    if(!is.numeric(x)){
      stop("Values in AdductMasses list must be numeric.")
    }
  }

  # It's weird if proton isn't included, so let the user know 
  if (any(c("Proton", "proton") %in% names(AdductMasses)) == FALSE) {
    warning("AdductMasses does not contain a mass for a proton. Proceeding with inputted masses.")
  }

  ###################
  ## LOAD GLOSSARY ##
  ###################

  # Load backend glossary
  Glossary <- data.table::fread(
    system.file("extdata", "Unimod_v20220602.csv", package = "pspecterlib")
  )

  ##################
  ## RUN ITERATOR ##
  ##################

  .calculate_molform_iterator <- function(Biomolecule, Name, Charge, AdductMass, AdductName) {

    #######################################
    ## GET MODIFICATION NAMES AND MASSES ##
    #######################################

    # If the input is a ProForma string
    if (BioType == "ProForma") {

      # Get modifications and sequence
      PTMs <- pspecterlib::convert_proforma(Biomolecule)

      # If PTMs is not a modifications_pspecter object (i.e. a sequence was just returned), use a divergent path
      if (!inherits(PTMs, "modifications_pspecter")) {

        # Get the molecular formula
        Formula <- pspecterlib::get_aa_molform(PTMs)

        # Format the molecular formula
        MolForm <- Formula %>% pspecterlib::collapse_molform()

        # There are no mass changes
        MassChanges <- 0

        # Get the monoisotopic mass
        MonoisotopicMass <- pspecterlib::get_monoisotopic(Formula)
        MonoMass <- (MonoisotopicMass + (Charge * AdductMass)) / Charge

        # Set MAI to NA if not included
        MAI <- NA

        # Add most abundant isotope
        if (AddMostAbundantIsotope) {

          # Get the isotope profile
          Isotopes <- pspecterlib::calculate_iso_profile(molform = Formula, algorithm = IsotopeAlgorithm, min_abundance = MinAbundance)

          # Get the most abundant isotope
          MAI <- Isotopes[which.max(Isotopes$abundance), "mass"] %>% unlist()
          MAI <- (MAI + (Charge * AdductMass)) / Charge

        }

        # Generate data table
        IsoMatchMS_MolForm <- data.table::data.table(
          "Biomolecules" = Biomolecule,
          "Identifiers" = Name,
          "Adduct Name" = AdductName,
          "Adduct Mass" = AdductMass,
          "Charge" = Charge,
          "Molecular Formula" = MolForm,
          "Mass Shift" = MassChanges,
          "Monoisotopic Mass" = MonoMass,
          "Most Abundant Isotope" = MAI
        )

        # Return object
        return(IsoMatchMS_MolForm)

      }

      # Pull values
      Sequence <- attr(PTMs, "pspecter")$cleaned_sequence
      Modifications <- attr(PTMs, "pspecter")$PTMs
      MassChanges <- attr(PTMs, "pspecter")$mass_changes

      ####################################
      ## Generate the Molecular Formula ##
      ####################################

      # Generate a pspecter molecule object
      Formula <- pspecterlib::get_aa_molform(Sequence)

      # If there are modifications, add those as well
      if (length(Modifications) > 0) {

        for (PTM in Modifications) {
          
          # Extract formula
          premolform <- Glossary[Glossary$Modification == PTM, 4:ncol(Glossary)] %>%
            dplyr::select(colnames(.)[!is.na(.)]) %>%
            paste0(colnames(.), ., collapse = "")

          # Add to formula
          Formula <- pspecterlib::add_molforms(Formula, pspecterlib::as.molform(premolform))
        }

      }

      #####################################################
      ## Add monoisotopic mass and most abundant isotope ##
      #####################################################

      # Calculate monoisotopic mass
      MonoisotopicMass <- pspecterlib::get_monoisotopic(Formula)

      # Get monoisotopic mass, which will always be the first one
      MonoMass <- (MonoisotopicMass + (Charge * AdductMass)) / Charge

      # Set MAI to NA if not included
      MAI <- NA

      # Add most abundant isotope
      if (AddMostAbundantIsotope) {

        # Get the isotope profile
        Isotopes <- pspecterlib::calculate_iso_profile(molform = Formula, algorithm = IsotopeAlgorithm, min_abundance = MinAbundance)

        # Get the most abundant isotope
        MAI <- Isotopes[which.max(Isotopes$abundance), "mass"] %>% unlist()
        MAI <- (MAI + (Charge * AdductMass)) / Charge

      }

      # Add mass changes if they exist
      if (length(MassChanges) > 0) {
        MonoMass <- MonoMass + (sum(MassChanges) / Charge)
        if (AddMostAbundantIsotope) {MAI <- MAI + (sum(MassChanges) / Charge)}
      }

      # Collapse the molecular formula
      MolForm <- Formula %>% pspecterlib::collapse_molform()

      # Generate data table
      IsoMatchMS_MolForm <- data.table::data.table(
        "Biomolecules" = Biomolecule,
        "Identifiers" = Name,
        "Adduct Name" = AdductName,
        "Adduct Mass" = AdductMass,
        "Charge" = Charge,
        "Molecular Formula" = MolForm,
        "Mass Shift" = sum(MassChanges),
        "Monoisotopic Mass" = MonoMass,
        "Most Abundant Isotope" = MAI
      )

      # Return object
      return(IsoMatchMS_MolForm)

    } else if (BioType == "Molecular Formula") {

      # Read the molecular formula
      Formula <- pspecterlib::as.molform(Biomolecule)

      # Get Monoisotopic Mass
      MonoisotopicMass <- pspecterlib::get_monoisotopic(Formula)

      # Adjust for the charge
      MonoMass <- (MonoisotopicMass + (Charge * AdductMass)) / Charge

      if (AddMostAbundantIsotope) {

        # Get the isotope profile
        Isotopes <- pspecterlib::calculate_iso_profile(molform = Formula, algorithm = IsotopeAlgorithm, min_abundance = MinAbundance)

        # Get the most abundant isotope
        MAI <- Isotopes[which.max(Isotopes$abundance), "mass"] %>% unlist()
        MAI <- (MAI + (Charge * AdductMass)) / Charge

      } else {MAI <- NA}

      IsoMatchMS_MolForm <- data.table::data.table(
        "Biomolecules" = Biomolecule,
        "Identifiers" = Name,
        "Adduct Name" = AdductName,
        "Adduct Mass" = AdductMass,
        "Charge" = Charge,
        "Molecular Formula" = Biomolecule,
        "Mass Shift" = 0,
        "Monoisotopic Mass" = MonoMass,
        "Most Abundant Isotope" = MAI

      )
    }

  }

  # Iterate through biomolecules and charges, and calculate all required values
  Pre_MolForms <- data.table::data.table(
    Pform = rep(Biomolecules, each = length(Charge)),
    theName = rep(Identifiers, each = length(Charge)),
    ChargeValue = Charge
  ) %>%
    unique()

  # Duplicating rows based in the number of AdductMasses
  Pre_MolForms <- Pre_MolForms[rep(seq_len(nrow(Pre_MolForms)), each = length(AdductMasses)),]

  # Adding the AdductMasses as a column, creating a row for every combo of Protein, Charge, and AdductMass
  Pre_MolForms$AdductMasses <- rep(as.numeric(AdductMasses), times = nrow(Pre_MolForms)/length(AdductMasses))
  
  # Add names
  Pre_MolForms$AdductNames <- rep(names(AdductMasses), times = nrow(Pre_MolForms)/length(AdductMasses))

  # Implement parallel computing for speed
  #doParallel::registerDoParallel(parallel::detectCores())
  #
  ## Collect results
  #All_MolForms <- foreach(it = 1:nrow(Pre_MolForms), .combine = "rbind") %dopar% {
  #  .calculate_molform_iterator(
  #    Biomolecule = Pre_MolForms$Pform[it],
  #    Name = Pre_MolForms$theName[it],
  #    Charge = Pre_MolForms$ChargeValue[it],
  #    AdductMass = Pre_MolForms$AdductMasses[it],
  #    AdductName = Pre_MolForms$AdductNames[it]
  #  )
  #}
  
  # Include an lapply option for debugging
  All_MolForms <- do.call(rbind, lapply(1:nrow(Pre_MolForms), function(it) {
    .calculate_molform_iterator(
      Biomolecule = Pre_MolForms$Pform[it],
      Name = Pre_MolForms$theName[it],
      Charge = Pre_MolForms$ChargeValue[it],
      AdductMass = Pre_MolForms$AdductMasses[it],
      AdductName = Pre_MolForms$AdductNames[it]
    )
  }))

  # Add class
  class(All_MolForms) <- c("IsoMatchMS_MolForm", class(All_MolForms))

  return(All_MolForms)

}
