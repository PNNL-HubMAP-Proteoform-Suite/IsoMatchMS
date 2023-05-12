context(test: match_biomolecule_to_ms1)

test_that("testing ms1 matching", {

  ### Download required test tiles ###

  InPeakData <- readRDS(system.file("testdata", "Intact_PeakData.RDS", package = "IsoMatchMS"))
  PepPeakData <- readRDS(system.file("testdata", "Peptides_PeakData.RDS", package = "IsoMatchMS"))
  PepMolforms <- readRDS(system.file("testdata", "Peptides_Molform.RDS", package = "IsoMatchMS"))
  InMolforms <- readRDS(system.file("testdata", "Intact_Molform.RDS", package = "IsoMatchMS"))

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
    "MolecularFormula must be a ProteoMatch_MolForm object."
  )

  expect_error(
    match_biomolecule_to_ms1(
      PeakData = PepPeakData,
      MolecularFormulas = PepMolforms,
      IsotopeRange = c(5, 20),
      PPMThreshold = 0,
      MinAbundance = 0.1,
    ),
    "PPMThreshold must be a nonzero numeric."
  )

  expect_error(
    match_biomolecule_to_ms1(
      PeakData = PepPeakData,
      MolecularFormulas = PepMolforms,
      IsotopeRange = c(5, 20),
      PPMThreshold = "wrong",
      MinAbundance = 0.1,
    ),
    "PPMThreshold must be a nonzero numeric."
  )

  expect_error(
    match_biomolecule_to_ms1(
      PeakData = PepPeakData,
      MolecularFormulas = PepMolforms,
      IsotopeRange = c(5, 20),
      PPMThreshold = 10,
      MinAbundance = 101,
    ),
    "MinAbundance should be a numeric between 0 and 100, inclusive."
  )

  expect_error(
    match_biomolecule_to_ms1(
      PeakData = PepPeakData,
      MolecularFormulas = PepMolforms,
      IsotopeRange = c(5, 20),
      PPMThreshold = 10,
      MinAbundance = "Wrong",
    ),
    "MinAbundance should be a numeric between 0 and 100, inclusive."
  )

  expect_error(
    match_biomolecule_to_ms1(
      PeakData = PepPeakData,
      MolecularFormulas = PepMolforms,
      IsotopeRange = c(5),
      PPMThreshold = 10,
      MinAbundance = 0.1,
    ),
    "IsotopeRange must be a numeric with at least two unique values."
  )

  # expect_error(
  #   match_biomolecule_to_ms1(
  #     PeakData = PepPeakData,
  #     MolecularFormulas = PepMolforms,
  #     MatchingAlgorithm = "closest peak",
  #     IsotopeRange = c("not numbers"),
  #     PPMThreshold = 10,
  #     MinAbundance = 0.1,
  #   ),
  #   "IsotopeRange must be a numeric with at least two unique values."
  # )

  ### Running Function ###

  #Peptide data with closest peak algorithm
  PepMatches <- match_biomolecule_to_ms1(
    PeakData = PepPeakData,
    MolecularFormulas = PepMolforms,
    IsotopeRange = c(5, 20),
    PPMThreshold = 10,
    MinAbundance = 0.1,
  )
  expect_true(inherits(PepMatches, "ProteoMatch_MatchedPeaks"))

  #Intact data with highest abundance algorithm
  InMatches <- match_biomolecule_to_ms1(
    PeakData = InPeakData,
    MolecularFormulas = InMolforms,
    IsotopeRange = c(5, 20),
    PPMThreshold = 10,
    MinAbundance = 0.1,
  )
  expect_true(inherits(InMatches, "ProteoMatch_MatchedPeaks"))



})
