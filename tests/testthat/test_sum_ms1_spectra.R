context(test: sum_ms1_spectra)

test_that("testing summing ms1 spectra", {

  ### Download required test tiles ###

  ### Testing inputs ###

  expect_error(
    sum_ms1_spectra(
      mzMLPath = system.file("testdata", "Intact_Protein_Summed_MS1.mzML", package = "ProteoMatch"),
      PPM_Round = "a",
      MinimumAbundance = 0.01,
      Messages = FALSE
      ),
    "PPM_Round must be a numeric of length 1"
  )

  expect_error(
    sum_ms1_spectra(
      mzMLPath = system.file("testdata", "Intact_Protein_Summed_MS1.mzML", package = "ProteoMatch"),
      PPM_Round = 0,
      MinimumAbundance = 0.01,
      Messages = FALSE
    ),
    "PPM_Round must be a numeric of length 1"
  )

  expect_error(
    sum_ms1_spectra(
      mzMLPath = system.file("testdata", "Intact_Protein_Summed_MS1.mzML", package = "ProteoMatch"),
      PPM_Round = 5,
      MinimumAbundance = 101,
      Messages = FALSE
    ),
    "MinimumAbundance should be a single numeric between 0 and 100."
  )

  expect_error(
    sum_ms1_spectra(
      mzMLPath = system.file("testdata", "Intact_Protein_Summed_MS1.mzML", package = "ProteoMatch"),
      PPM_Round = 5,
      MinimumAbundance = 0.01,
      Messages = "Wrong"
    ),
    "Messages must be TRUE or FALSE"
  )

  ### Running Function ###

  summedSpec <- sum_ms1_spectra(
    mzMLPath = system.file("testdata", "20220228_VU_pancreas_40um_IMAGE-16078AVG.mzML", package = "ProteoMatch"),
    PPM_Round = 5,
    MinimumAbundance = 0.01,
    Messages = TRUE
  )
  expect_true(inherits(summedSpec, "peak_data"))

})
