context(test: run_proteomatch)

test_that("testing proteomatch wrapper function", {

  ### Download required test tiles ###

  PepPeakData <- readRDS(system.file("testdata", "Intact_PeakData.RDS", package = "ProteoMatch"))
  InPeakData <- readRDS(system.file("testdata", "Peptides_PeakData.RDS", package = "ProteoMatch"))
  InPeakList <- read.csv(system.file("testdata", "Intact_Proteins_List_Short.csv", package = "ProteoMatch"))
  PepPeakList <- read.csv(system.file("testdata", "Peptides_List_Short.csv", package = "ProteoMatch"))

  InMols <- readRDS(system.file("testdata", "Intact_Molform.RDS", package = "ProteoMatch"))

  ### Testing inputs ###

  expect_error(
    run_proteomatch(
      Modifications = InPeakList$Formula,
      ModType = "ProForma",
      SummedSpectra = InPeakData,
      SettingsFile = system.file("testdata", "Peptides_List_Short.csv", package = "ProteoMatch"),
      Proteins = InPeakList$Protein
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
    run_proteomatch(
      Modifications = InPeakList$Formula,
      ModType = "ProForma",
      SummedSpectra = InPeakData,
      SettingsFile = system.file("testdata", "Intact_Protein_Defaults.xlsx", package = "ProteoMatch"),
      Proteins = InPeakList$Protein,
      Messages = "Wrong"
    ),
    "Messages should be a TRUE or FALSE."
  )

  ### Running Function ###

  # trel1 <- run_proteomatch(
  #   Modifications = InPeakList$Proteoform,
  #   ModType = "ProForma",
  #   SummedSpectra = InPeakData,
  #   SettingsFile = system.file("testdata", "Intact_Protein_Defaults.xlsx", package = "ProteoMatch"),
  #   Proteins = InPeakList$Protein,
  #   Messages = F
  # )
  #
  # trel2 <- run_proteomatch(
  #   Modifications = InMols$`Molecular Formula`,
  #   ModType = "Molecular Formula",
  #   SummedSpectra = InPeakData,
  #   SettingsFile = system.file("testdata", "Intact_Protein_Defaults.xlsx", package = "ProteoMatch"),
  #   Proteins = NULL,
  #   Messages = TRUE
  # )

})
