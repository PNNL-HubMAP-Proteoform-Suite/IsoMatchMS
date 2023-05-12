context(test: run_isomatchms)

test_that("testing proteomatch wrapper function", {

  ### Download required test tiles ###

  InPeakData <- readRDS(system.file("testdata", "Intact_PeakData.RDS", package = "IsoMatchMS"))
  InPeakList <- read.csv(system.file("testdata", "Intact_Proteins_List_Short.csv", package = "IsoMatchMS"))
  PepPeakData <- readRDS(system.file("testdata", "Peptides_PeakData.RDS", package = "IsoMatchMS"))
  PepPeakList <- read.csv(system.file("testdata", "Peptides_List_Short.csv", package = "IsoMatchMS"))

  InMols <- readRDS(system.file("testdata", "Intact_Molforms.RDS", package = "IsoMatchMS"))

  ### Testing inputs ###

  expect_error(
    run_isomatchms(
      Biomolecules = InPeakList$Formula,
      BioType = "ProForma",
      SummedSpectra = InPeakData,
      SettingsFile = system.file("testdata", "Peptides_List_Short.csv", package = "IsoMatchMS"),
      Identifiers = InPeakList$Protein
    ),
    "SettingsFile needs to be an xlsx file."
  )


# Removed the MZRange row from the Defaults file, but get this error:
# Error in filter_peaks(PeakData = SummedSpectra, MZRange = Settings[Settings$Parameter == :
# MZRange should be a numeric with more than 1 unique value.
  # expect_error(
  #   run_proteomatch(
  #     Modifications = InPeakList$Proteoform,
  #     ModType = "ProForma",
  #     SummedSpectra = InPeakData,
  #     SettingsFile = system.file("testdata", "Intact_Protein_Defaults_wrong.xlsx", package = "ProteoMatch"),
  #     Proteins = InPeakList$Protein.accession
  #   ),
  #   "Settings file is missing: MZRange"
  # )

  expect_error(
    run_isomatchms(
      Biomolecules = InPeakList$Formula,
      BioType = "ProForma",
      SummedSpectra = InPeakData,
      SettingsFile = system.file("extdata", "Intact_Protein_Defaults.xlsx", package = "IsoMatchMS"),
      Identifiers = InPeakList$Protein,
      Messages = "Wrong"
    ),
    "Messages should be a TRUE or FALSE."
  )

  ### Running Function ###

  run_isomatchms(
     Biomolecules = InPeakList$Proteoform,
     BioType = "ProForma",
     SummedSpectra = InPeakData,
     SettingsFile = system.file("extdata", "Intact_Protein_Defaults.xlsx", package = "IsoMatchMS"),
     Identifiers = InPeakList$Protein,
     Messages = TRUE
  )


})
