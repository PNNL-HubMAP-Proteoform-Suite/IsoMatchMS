context("test: plot_Ms1Match")

test_that("testing plotting of IsoMatchMS_MatchedPeaks object", {

  ### Download required test tiles ###

  PepPeakData <- readRDS(system.file("testdata", "Intact_PeakData.RDS", package = "IsoMatchMS"))
  InPeakData <- readRDS(system.file("testdata", "Peptides_PeakData.RDS", package = "IsoMatchMS"))
  PepMatches <- readRDS(system.file("testdata", "Peptide_Matches.RDS", package = "IsoMatchMS"))
  InMatches <- readRDS(system.file("testdata", "Intact_Matches.RDS", package = "IsoMatchMS"))


  ### Testing inputs ###

  expect_error(
    plot_Ms1Match(
      PeakData = "Wrong",
      Ms1Match = PepMatches,
      ID = 1,
      Window = 5
    ),
    "PeakData must be a pspecterlib peak_data object."
  )

  expect_error(
    plot_Ms1Match(
      PeakData = PepPeakData,
      Ms1Match = "Wrong",
      ID = 1,
      Window = 5
    ),
    "Ms1Match must be a IsoMatchMS_MatchedPeaks object from match_biomolecule_to_ms1"
  )

  expect_error(
    plot_Ms1Match(
      PeakData = PepPeakData,
      Ms1Match = PepMatches,
      ID = c(1, 2),
      Window = 5
    ),
    "ID should be of length 1."
  )

  expect_error(
    plot_Ms1Match(
      PeakData = PepPeakData,
      Ms1Match = PepMatches,
      ID = 100,
      Window = 5
    ),
    "100 is not a recognized ID."
  )

  expect_error(
    plot_Ms1Match(
      PeakData = PepPeakData,
      Ms1Match = PepMatches,
      ID = 1,
      Window = "a"
    ),
    "Window must be a single numeric."
  )

  expect_error(
    plot_Ms1Match(
      PeakData = PepPeakData,
      Ms1Match = PepMatches,
      ID = 1,
      Window = c(5, 6)
    ),
    "Window must be a single numeric."
  )

  ### Running Function ###

  #Peptide data
  # Run two examples with two charge states
  MolForms_Test <- calculate_molform(
    Biomolecule = c("M.SS[Methyl]S.V", "M.SS[6]S[7].V"),
    BioType = "ProForma",
    Identifiers = c("Test1", "Test2"),
    Charge = 1:2
  )
  
  # Generate some experimental peak data to match
  PeakData <- pspecterlib::make_peak_data(
    MZ = c(147.5684, 148.0699, 148.5708, 149.0721, 149.5731,
           294.1296, 295.1325, 296.1343, 297.1369, 298.1390),
    Intensity = c(868.3680036, 110.9431876, 18.7179196, 1.7871629, 0.1701294,
                  868.3680036, 110.9431876, 18.7179196, 1.7871629, 0.1701294)
  )
  
  # Run algorithm
  IsoMatch <- match_biomolecule_to_ms1(
    PeakData = PeakData,
    MolecularFormula = MolForms_Test,
    IsotopeMinimum = 3
  )
  
  # Make plot
  InPlot <- plot_Ms1Match(
    PeakData = PeakData,
    Ms1Match = IsoMatch,
    ID = 1
  )

  expect_true(inherits(InPlot, "plotly"))



})
