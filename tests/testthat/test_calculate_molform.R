context("test: calculate_molform")

test_that("testing molform calculation", {

  ### Download required test tiles ###

  Peptide_data <- read.csv(system.file("testdata", "Peptides_List_Short.csv", package = "IsoMatchMS"))
  Protein_data <- read.csv(system.file("testdata", "Intact_Proteins_List_Short.csv", package = "IsoMatchMS"))
  InMols <- readRDS(system.file("testdata", "Intact_Molforms.RDS", package = "IsoMatchMS"))

  ### Testing inputs ###

  # Biomolecules should be strings
  expect_error(
    calculate_molform(Biomolecules = 4),
    "Biomolecules must be a vector of characters."
  )

  # Blank or empty biomolecules are not permitted
  expect_error(
    calculate_molform(Biomolecules = c("AVEEMITITK", NA, "DQVANSAFVER"), BioType = "ProForma"),
    "Biomolecules cannot be NULL, NA, or blank."
  )

  # Check that BioType is "Molecular Formula" or "ProForma"
  expect_error(
    calculate_molform(Biomolecules = c("AVEEMITITK", "DQVANSAFVER"), BioType = "wrong"),
    "BioType can either be 'ProForma' or 'Molecular Formula'."
  )

  # Identifiers should be a string
  expect_error(
    calculate_molform(Biomolecules = c("AVEEMITITK", "DQVANSAFVER"), BioType = "ProForma", Identifiers = 4),
    "Identifiers must be a vector of characters."
  )

  # Proteoform and Identifiers should be the same length
  expect_error(
    calculate_molform(Biomolecules = c("AVEEMITITK", "DQVANSAFVER"), BioType = "ProForma", Identifiers = c("P12839")),
    "Biomolecules and Identifiers must be of the same length."
  )

  # Charge must be numeric and will be rounded to integers
  expect_error(
    calculate_molform(Biomolecules = c("AVEEMITITK", "DQVANSAFVER"), BioType = "ProForma", Charge = "wrong"),
    "Charge should be numeric."
  )

  # Check that AddMostAbundantIsotope is a TRUE/FALSE
  expect_error(
    calculate_molform(Biomolecules = c("AVEEMITITK", "DQVANSAFVER"), BioType = "ProForma", AddMostAbundantIsotope = "wrong"),
    "AddMostAbundnantIsotope must be true or false."
  )

  # Adduct mass must be a single numeric
  expect_error(
    calculate_molform(Biomolecules = c("AVEEMITITK", "DQVANSAFVER"), BioType = "ProForma", AdductMasses = "wrong"),
    "Values in AdductMasses list must be numeric."
  )

  # Adduct mass cannot be longer than 5 
  expect_error(
    calculate_molform(Biomolecules = c("AVEEMITITK", "DQVANSAFVER"), BioType = "ProForma", 
                      AdductMasses = c("Test1" = 1, "Test2" = 2, "Test3" = 3, "Test4" = 4, "Test5" = 5, "Test6" = 6)),
    "AdductMasses must be a list of length 1 to 5."
  )
  
  # Warning if no proton is included
  expect_warning(
    calculate_molform(Biomolecules = c("AVEEMITITK", "DQVANSAFVER"), BioType = "ProForma", 
                      AdductMasses = c("Test1" = 1)),
    "AdductMasses does not contain a mass for a proton. Proceeding with inputted masses."
  )
  
  
  ### Function Tests ###

  # Testing peptide data
  Pep_mols <- calculate_molform(Biomolecules = Peptide_data$Proteoform, BioType = "ProForma", Identifiers = Peptide_data$Protein, AddMostAbundantIsotope = TRUE)
  expect_true(inherits(Pep_mols, "IsoMatchMS_MolForm"))
  
  # Testing protein data
  Prot_mols1 <- calculate_molform(Biomolecules = Protein_data$Proteoform, BioType = "ProForma", Identifiers = Protein_data$Protein.accession, AddMostAbundantIsotope = TRUE)
  expect_true(inherits(Prot_mols1, "IsoMatchMS_MolForm"))
  Prot_mols2 <- calculate_molform(Biomolecules = InMols$`Molecular Formula`, BioType = "Molecular Formula", Identifiers = Protein_data$Protein, AddMostAbundantIsotope = TRUE)
  expect_true(inherits(Prot_mols2, "IsoMatchMS_MolForm"))
  Prot_mols3 <- calculate_molform(Biomolecules = InMols$`Molecular Formula`, BioType = "Molecular Formula", Identifiers = Protein_data$Protein, AddMostAbundantIsotope = FALSE)
  expect_true(inherits(Prot_mols3, "IsoMatchMS_MolForm"))
  
})
