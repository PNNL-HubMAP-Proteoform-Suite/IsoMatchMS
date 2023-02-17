#' Generate a Molecular Formula Table from provided molecular formulas, or ProForma strings
#' 
#' @description Returns the "IsoMatchMS_MolForm" object with each biomolecule's molecular formula and adjusted mass.
#'
#' @param Biomolecules A vector of Molecular Formulas or ProForma strings (i.e. "M.AA`[`Acetyl`]`AA`[`3.2`]`.V"). 
#'    ProForma strings can be pulled from mzid files (MS-GF+, MSPathFinder) 
#'    with pull_modifications_from_mzid, or made with create_proforma for MSPathFinder,
#'    ProSight, and pTop modifications. TopPIC proteoforms are provided as ProForma
#'    strings. Required.
#' @param BioType A string indicating whether the Modifications are "ProForma" strings or "Molecular Formula". Required. 
#' @param Identifiers A vector of identifiers for each biomolecule (i.e. protein, glycan, etc.). Optional.
#' @param Charge The range of charges to test. Default is 1.
#' @param AddMostAbundantIsotope A flag to determine whether the Most Abundant Isotope (MAI) should be calculated for 
#'     every function. This parameter will slow down tool. Default is FALSE. 
#' @param MinAbundance The minimum abundance (calculated intensity) permitted to be matched. 
#'     Default is 0.1, which is 0.1%. Used for most abundant isotope. This is a pspecterlib-specific
#'     parameter and shouldn't need to be changed for IsoMatchMS. 
#' @param AdductMasses The masses of adducts to be tested. Proton Adducts are the default. 
#'     Default is c(1.00727647).
#'
#' @details
#' The data.table outputted by this function returns 8 columns.
#' \tabular{ll}{
#' Biomolecules \tab The provided biomolecule molecular formulas or proforma strings \cr
#' \tab \cr
#' Identifiers \tab The provided biomolecule identifier \cr
#' \tab \cr
#' Adduct \tab The provided adduct \cr
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
#' @returns A IsoMatchMS_MolForm object, which is a data.table containing the
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
#'    Charge = 1:3
#' )
#'
#' # Run two examples with two charge states
#' calculate_molform(
#'    Biomolecules = c("M.SS[Methyl]S.V", "M.S[Methyl]S[22]S[23].V"),
#'    BioType = "ProForma",
#'    Charge = 1:2
#' )
#'
#' # Run an example with molecular formulas
#' calculate_molform(
#'    Biomolecules = c("C6H12O6", "C2H4O1"),
#'    BioType = "Molecular Formula",
#'    Identifiers = c("Glucose", "Acetyl")
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
                              MinAbundance = 0.1,
                              AdductMasses = c(1.00727647)) {

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
    stop("BioType can either be 'ProForma' or 'Molecular Formula'")
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
  
  # Proton Mass must be a single numeric
  if (length(AdductMasses) != 1 | !is.numeric(AdductMasses)) {
    stop("AdductMasses must be a single numeric.")
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

  .calculate_molform_iterator <- function(Biomolecule, Name, Charge) {

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
        MolForm <- Formula %>% 
          .[. > 0] %>% 
          paste0(names(.), ., collapse = "")
        
        # There are no mass changes
        MassChanges <- 0
        
        # Get the monoisotopic mass
        MonoisotopicMass <- pspecterlib::get_monoisotopic(Formula)
        MonoMass <- (MonoisotopicMass + (Charge * AdductMasses)) / Charge
        
        # Set MAI to NA if not included 
        MAI <- NA
        
        # Add most abundant isotope
        if (AddMostAbundantIsotope) {
          
          # Get the isotope profile 
          Isotopes <- pspecterlib::calculate_iso_profile(Formula, MinAbundance)
          
          # Get the most abundant isotope
          MAI <- Isotopes[which.max(Isotopes$abundance), "mass"] %>% unlist()
          MAI <- (MAI + (Charge * AdductMasses)) / Charge
          
        }
        
      
        # Generate data table
        ProteoMatch_MolForm <- data.table::data.table(
          "Biomolecules" = Biomolecule,
          "Name" = Name, 
          "Charge" = Charge,
          "Molecular Formula" = MolForm, 
          "Mass Shift" = MassChanges, 
          "Monoisotopic Mass" = MonoMass,
          "Most Abundant Isotope" = MAI
        )
        
        # Return object
        return(ProteoMatch_MolForm)
        
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
          premolform <- Glossary[Glossary$Modification == PTM, c(4:ncol(Glossary))] %>%
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
      MonoMass <- (MonoisotopicMass + (Charge * AdductMasses)) / Charge
      
      # Set MAI to NA if not included 
      MAI <- NA
      
      # Add most abundant isotope
      if (AddMostAbundantIsotope) {
        
        # Get the isotope profile 
        Isotopes <- pspecterlib::calculate_iso_profile(Formula, MinAbundance)
        
        # Get the most abundant isotope
        MAI <- Isotopes[which.max(Isotopes$abundance), "mass"] %>% unlist()
        MAI <- (MAI + (Charge * AdductMasses)) / Charge
        
      }
  
      # Add mass changes if they exist
      if (length(MassChanges) > 0) {
        MonoMass <- MonoMass + (sum(MassChanges) / Charge)
        if (AddMostAbundantIsotope) {MAI <- MAI + (sum(MassChanges) / Charge)}
      }
      
      # Collapse the molecular formula 
      CleanedNames <- Formula[Formula > 0]
      MolForm <- paste0(names(CleanedNames), CleanedNames, collapse = "")
      
      # Generate data table
      ProteoMatch_MolForm <- data.table::data.table(
        "Biomolecules" = Biomolecule,
        "Identifiers" = Name, 
        "Adduct" = AdductMasses,
        "Charge" = Charge,
        "Molecular Formula" = MolForm, 
        "Mass Shift" = sum(MassChanges), 
        "Monoisotopic Mass" = MonoMass,
        "Most Abundant Isotope" = MAI
      )
  
      # Return object
      return(ProteoMatch_MolForm)
  
    } else if (BioType == "Molecular Formula") {
      
      # Read the molecular formula
      Formula <- pspecterlib::as.molform(Biomolecule)
      
      # Get Monoisotopic Mass
      MonoisotopicMass <- pspecterlib::get_monoisotopic(Formula)
      
      # Adjust for the charge
      MonoMass <- (MonoisotopicMass + (Charge * AdductMasses)) / Charge
      
      if (AddMostAbundantIsotope) {
      
        # Get the isotope profile 
        Isotopes <- pspecterlib::calculate_iso_profile(Formula, MinAbundance)
        
        # Get the most abundant isotope
        MAI <- Isotopes[which.max(Isotopes$abundance), "mass"] %>% unlist()
        MAI <- (MAI + (Charge * AdductMasses)) / Charge
        
      } else {MAI <- NA}
      
      ProteoMatch_MolForm <- data.table::data.table(
        "Biomolecules" = Biomolecules,
        "Identifiers" = Name,
        "Adduct" = AdductMasses,
        "Charge" = Charge,
        "Molecular Formula" = Biomolecules, 
        "Mass Shift" = 0,
        "Monoisotopic Mass" = MonoMass,
        "Most Abundant Isotope" = MAI
        
      )
    }
    
  }

  # Iterate through Biomoleculess and charges, and calculate all required values
  Pre_MolForms <- data.table::data.table(
    Pform = rep(Biomolecules, each = length(Charge)),
    theName = rep(Identifiers, each = length(Charge)),
    ChargeValue = Charge
  ) %>%
    unique()
  
  # Implement parallel computing for speed 
  doParallel::registerDoParallel(parallel::detectCores())
  
  # Collect results
  All_MolForms <- foreach(it = 1:nrow(Pre_MolForms), .combine = rbind) %dopar% {
    .calculate_molform_iterator(
      Biomolecule = Pre_MolForms$Pform[it],
      Name = Pre_MolForms$theName[it],
      Charge = Pre_MolForms$ChargeValue[it]
    )
  }

  # Add class
  class(All_MolForms) <- c("IsoMatchMS_MolForm", class(All_MolForms))

  return(All_MolForms)

}
