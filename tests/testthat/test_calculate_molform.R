context(test: calculate_molform)

test_that("testing molform calculation", {

  ### Download required test tiles ###

  Peptide_data <- read.csv(system.file("testdata", "Peptides_List_Short.csv", package = "ProteoMatch"))
  Protein_data <- read.csv(system.file("testdata", "Intact_Proteins_List_Short.csv", package = "ProteoMatch"))
  InMols <- readRDS(system.file("testdata", "Intact_Molform.RDS", package = "ProteoMatch"))

  ### Testing inputs ###

  # Modifications should be a string
  expect_error(
    calculate_molform(Modifications = 4, ModType = "ProForma"),
    "Modifications must be a vector of characters."
  )

  # Blank or empty modifications are not permitted
  expect_error(
    calculate_molform(Modifications = c("AVEEMITITK", NA, "DQVANSAFVER"), ModType = "ProForma"),
    "Modifications cannot be NULL, NA, or blank."
  )

  # Check that ModType is "Molecular Formula" or "ProForma"
  expect_error(
    calculate_molform(Modifications = c("AVEEMITITK", "DQVANSAFVER"), ModType = "wrong"),
    "ModType can either be 'ProForma' or 'Molecular Formula'"
  )

  # Proteins should be a string
  expect_error(
    calculate_molform(Modifications = c("AVEEMITITK", "DQVANSAFVER"), ModType = "ProForma", Proteins = 4),
    "Proteins must be a vector of characters."
  )

  # Proteoform and Proteins should be the same length
  expect_error(
    calculate_molform(Modifications = c("AVEEMITITK", "DQVANSAFVER"), ModType = "ProForma", Proteins = c("P12839")),
    "Modifications and Proteins must be of the same length."
  )

  # Charge must be numeric and will be rounded to integers
  expect_error(
    calculate_molform(Modifications = c("AVEEMITITK", "DQVANSAFVER"), ModType = "ProForma", Charge = "wrong"),
    "Charge should be numeric."
  )

  # Check that AddMostAbundantIsotope is a TRUE/FALSE
  expect_error(
    calculate_molform(Modifications = c("AVEEMITITK", "DQVANSAFVER"), ModType = "ProForma", AddMostAbundantIsotope = "wrong"),
    "AddMostAbundnantIsotope must be true or false."
  )

  # Proton Mass must be a single numeric
  expect_error(
    calculate_molform(Modifications = c("AVEEMITITK", "DQVANSAFVER"), ModType = "ProForma", ProtonMass = "wrong"),
    "ProtonMass must be a single numeric."
  )

  ### Running Function ###

  #testing peptide data
  Pep_mols <- calculate_molform(Modifications = Peptide_data$Proteoform, ModType = "ProForma", Proteins = Peptide_data$Protein)
  expect_true(inherits(Pep_mols, "ProteoMatch_MolForm"))

  #testing intact protein data
  # Prot_mols <- calculate_molform(Modifications = Protein_data$Proteoform, ModType = "ProForma", Proteins = Protein_data$Protein, AddMostAbundantIsotope = F)
  # expect_true(inherits(Prot_mols, "ProteoMatch_MolForm"))

  Prot_mols2 <- calculate_molform(Modifications = InMols$`Molecular Formula`, ModType = "Molecular Formula", Proteins = Protein_data$Protein, AddMostAbundantIsotope = F)
  expect_true(inherits(Prot_mols2, "ProteoMatch_MolForm"))
})
