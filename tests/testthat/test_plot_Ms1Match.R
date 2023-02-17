context(test: plot_Ms1Match)

test_that("testing plotting of ProteoMatch_MatchedPeaks object", {

  ### Download required test tiles ###

  PepPeakData <- readRDS(system.file("testdata", "Intact_PeakData.RDS", package = "ProteoMatch"))
  InPeakData <- readRDS(system.file("testdata", "Peptides_PeakData.RDS", package = "ProteoMatch"))
  PepMatches <- readRDS(system.file("testdata", "Peptide_Matches.RDS", package = "ProteoMatch"))
  InMatches <- readRDS(system.file("testdata", "Intact_Matches.RDS", package = "ProteoMatch"))


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
    "Ms1Match must be a ProteoMatch_MatchedPeaks object from match_proteoform_to_ms1"
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
  PepPlot <- plot_Ms1Match(
    PeakData = PepPeakData,
    Ms1Match = PepMatches,
    ID = 1,
    Window = 5
  )
  expect_true(inherits(PepPlot, "ggplot"))

  #Intact data
  InPlot <- plot_Ms1Match(
    PeakData = InPeakData,
    Ms1Match = InMatches,
    ID = 2,
    Window = 5
  )
  expect_true(inherits(InPlot, "ggplot"))



})
