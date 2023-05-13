context("test: match_biomolecule_to_ms1")

test_that("testing ms1 matching", {

  ### Download required test tiles ###

  InPeakData <- readRDS(system.file("testdata", "Intact_PeakData.RDS", package = "IsoMatchMS"))
  PepPeakData <- readRDS(system.file("testdata", "Peptides_PeakData.RDS", package = "IsoMatchMS"))
  PepMolforms <- readRDS(system.file("testdata", "Peptides_Molforms.RDS", package = "IsoMatchMS"))
  InMolforms <- readRDS(system.file("testdata", "Intact_Molforms.RDS", package = "IsoMatchMS"))

  ### Testing inputs ###

  expect_error(
    match_biomolecule_to_ms1(
      PeakData = "Wrong",
      MolecularFormulas = PepMolforms
    ),
    "PeakData must be a pspecterlib peak_data object."
  )

  expect_error(
    match_biomolecule_to_ms1(
      PeakData = PepPeakData,
      MolecularFormulas = "Wrong"
    ),
    "MolecularFormula must be an IsoMatchMS_MolForm object."
  )

  expect_error(
    match_biomolecule_to_ms1(
      PeakData = PepPeakData,
      MolecularFormulas = PepMolforms,
      AbundanceThreshold = "Wrong"
    ),
    "AbundanceThreshold must be a nonzero numeric."
  )

  expect_error(
    match_biomolecule_to_ms1(
      PeakData = PepPeakData,
      MolecularFormulas = PepMolforms,
      IsotopeMinimum = 3,
      PPMThreshold = 10,
      IsotopeAlgorithm = "Wrong"
    ),
    "IsotopeAlgorithm must either be 'isopat' or 'Rdisop'."
  )
  
  expect_error(
    match_biomolecule_to_ms1(
      PeakData = PepPeakData,
      MolecularFormulas = PepMolforms,
      PPMThreshold = "Wrong"
    ),
    "PPMThreshold must be a nonzero numeric."
  )
  
  expect_error(
    match_biomolecule_to_ms1(
      PeakData = PepPeakData,
      MolecularFormulas = PepMolforms,
      IsotopeMinimum = "Wrong"
    ),
    "IsotopeRange must be a non-zero numeric value."
  )
  
  expect_error(
    match_biomolecule_to_ms1(
      PeakData = PepPeakData,
      MolecularFormulas = PepMolforms,
      IsotopeMinimum = 20
    ),
    "No fragmentation patterns found. Try increasing the PPMThreshold, decreasing the minimum isotope range, lowering the noise filter, or lowering the minimum abundance."
  )

  ### Running Function ###

  # Peptide data with isopat isotoping algorithm
  PepMatches <- match_biomolecule_to_ms1(
    PeakData = PepPeakData,
    MolecularFormulas = PepMolforms,
    IsotopeAlgorithm = "isopat"
  )
  expect_true(inherits(PepMatches, "IsoMatchMS_MatchedPeaks"))

  # Intact data with Rdisop isotoping algorithm
  InMatches <- match_biomolecule_to_ms1(
    PeakData = InPeakData,
    MolecularFormulas = InMolforms
  )
  expect_true(inherits(InMatches, "IsoMatchMS_MatchedPeaks"))



})
