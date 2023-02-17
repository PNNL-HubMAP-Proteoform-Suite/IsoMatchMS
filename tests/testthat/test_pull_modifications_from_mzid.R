context(test: pull_modifications_from_mzid)

test_that("testing pulling of modification data from mzid", {

  ### Download required test tiles ###

  # Create a temporary directory and copy example data there
  tmpdir <- tempdir()

  file <- "https://raw.githubusercontent.com/EMSL-Computing/PSpecteR/master/pspecter_container/TestFiles/BottomUp/BottomUp.mzid"

  download.file(file, file.path(tmpdir, tail(unlist(strsplit(file, "/")), 1)))

  ### Testing inputs ###

  expect_error(
    pull_modifications_from_mzid("/This/Doesnt/Exist.mzID"),
    "ID file path must exist."
  )

  expect_error(
    pull_modifications_from_mzid("./DESCRIPTION"),
    "ID file must be an mzid file."
  )

  ### Running Function ###

  #testing BottomUp example
  BotMods <- pull_modifications_from_mzid(file.path(tmpdir, "BottomUp.mzid"))
  expect_true("Proteoform" %in% colnames(BotMods))
})
