context(test: filter_peak)

test_that("testing peak filtering", {

  ### Download required test tiles ###

  PepPeakData <- readRDS(system.file("testdata", "Peptides_PeakData.RDS", package = "IsoMatchMS"))
  InPeakData <- readRDS(system.file("testdata", "Intact_PeakData.RDS", package = "IsoMatchMS"))

  ### Testing inputs ###

  expect_error(
    filter_peaks(
      PeakData = "Wrong",
      MZRange = c(500, 1550),
      NoiseFilter = 2
    ),
    "PeakData should be a peak_data object from the pspecterlib get_peak_data or make_peak_data functions."
  )

  expect_message(
    filter_peaks(
      PeakData = PepPeakData,
      MZRange = c(500, 1550, 1200),
      NoiseFilter = 2
    ),
    "More than 2 values found in MZRange. Only the min and max will be used."
  )
  
  
  expect_error(
    filter_peaks(
      PeakData = PepPeakData,
      MZRange = "Wrong",
      NoiseFilter = 2
    ),
    "MZRange should be a numeric with more than 1 unique value."
  )

  expect_error(
    filter_peaks(
      PeakData = PepPeakData,
      MZRange = c(500, 1550),
      NoiseFilter = "Wrong"
    ),
    "NoiseFilter should be a numeric between 0 and 100, inclusive."
  )

  ### Running Function ###

  FiltPep <- filter_peaks(
    PeakData = PepPeakData,
    MZRange = c(0, 1800),
    NoiseFilter = 1
  )
  expect_true(inherits(FiltPep, "peak_data"))

  FiltIn <- filter_peaks(
    PeakData = InPeakData,
    MZRange = c(5000, 6000),
    NoiseFilter = 2
  )
  expect_true(inherits(FiltIn, "peak_data"))


})
