context("test: create_proforma")

test_that("testing proforma string creation", {

  ### Testing inputs ###

  expect_error(
    create_proforma(
      Sequence = 4,
      Modifications = c("Methyl 8,Phospho 12,Phospho 22,DiMethyl 30", "TriMethyl 6,Phospho 15", ""),
      Tool = "MSPathFinder",
    ),
    "Sequence should be a vector of character strings."
  )

  expect_error(
    create_proforma(
      Sequence = c("Th1s 1sn'7 a S3qu3nc3"),
      Modifications = c("Methyl 8,Phospho 12,Phospho 22,DiMethyl 30", "TriMethyl 6,Phospho 15", ""),
      Tool = "MSPathFinder",
    ),
    "Not all sequences in Sequence are valid. The first example is: Th1s 1sn'7 a S3qu3nc3"
  )

  expect_error(
    create_proforma(
      Sequence = c("GRGKTGGKARAKAKSRSSRAGLQFPVGRVHRL", "KKTRIIPRHLQLAIRNDEELNKLLGGVTIAY", "TEST"),
      Modifications = 4,
      Tool = "MSPathFinder",
    ),
    "Modifications must be a vector of character strings."
  )

  expect_error(
    create_proforma(
      Sequence = c("GRGKTGGKARAKAKSRSSRAGLQFPVGRVHRL", "KKTRIIPRHLQLAIRNDEELNKLLGGVTIAY", "TEST"),
      Modifications = c("Methyl 8,Phospho 12,Phospho 22,DiMethyl 30"),
      Tool = "MSPathFinder",
    ),
    "Sequence and Modifications should be the same length."
  )

  expect_error(
    create_proforma(
      Sequence = c("GRGKTGGKARAKAKSRSSRAGLQFPVGRVHRL", "KKTRIIPRHLQLAIRNDEELNKLLGGVTIAY", "TEST"),
      Modifications = c("Methyl 8,Phospho 12,Phospho 22,DiMethyl 30", "TriMethyl 6,Phospho 15", ""),
      Tool = "Hammer",
    ),
    "Hammer is not a recognized type of Tool. If you would like to add it, please make an issue on the github page."
  )

  expect_error(
    create_proforma(
      Sequence = c("(49)M(37)SGRGKQG", "SG(67)RGKQGGKARAKAKSRSSRAG", "TEST"),
      Modifications = c("N-acetyl-L-methionine (49), O-phospho-L-serine (37)", "omega-N,omega-N'-dimethyl-L-arginine (67)", ""),
      Tool = "ProSight",
    ),
    "ConversionList is required for ProSight datasets."
  )

  expect_error(
    create_proforma(
      Sequence = c("GRGKTGGKARAKAKSRSSRAGLQFPVGRVHRL", "KKTRIIPRHLQLAIRNDEELNKLLGGVTIAY", "TEST"),
      Modifications = c("Methyl 8,Phospho 12,Phospho 22,DiMethyl 30", "TriMethyl 6,Phospho 15", ""),
      Tool = "MSPathFinder",
      Scan = c("not", "right", "Scan"),
    ),
    "Scan should be a vector of numeric scan numbers."
  )

  expect_error(
    create_proforma(
      Sequence = c("GRGKTGGKARAKAKSRSSRAGLQFPVGRVHRL", "KKTRIIPRHLQLAIRNDEELNKLLGGVTIAY", "TEST"),
      Modifications = c("Methyl 8,Phospho 12,Phospho 22,DiMethyl 30", "TriMethyl 6,Phospho 15", ""),
      Tool = "MSPathFinder",
      Scan = c(1, 2),
    ),
    "Sequence and Scan should be the same length."
  )

  expect_error(
    create_proforma(
      Sequence = c("GRGKTGGKARAKAKSRSSRAGLQFPVGRVHRL", "KKTRIIPRHLQLAIRNDEELNKLLGGVTIAY", "TEST"),
      Modifications = c("Methyl 8,Phospho 12,Phospho 22,DiMethyl 30", "TriMethyl 6,Phospho 15", ""),
      Tool = "MSPathFinder",
      Protein = 4
    ),
    "Protein should be a vector of characters."
  )

  expect_error(
    create_proforma(
      Sequence = c("GRGKTGGKARAKAKSRSSRAGLQFPVGRVHRL", "KKTRIIPRHLQLAIRNDEELNKLLGGVTIAY", "TEST"),
      Modifications = c("Methyl 8,Phospho 12,Phospho 22,DiMethyl 30", "TriMethyl 6,Phospho 15", ""),
      Tool = "MSPathFinder",
      Protein = c("Protein33")
    ),
    "Sequence and Protein should be the same length."
  )

  expect_error(
    create_proforma(
      Sequence = c("(49)M(37)SGRGKQG", "SG(67)RGKQGGKARAKAKSRSSRAG", "TEST"),
      Modifications = c("N-acetyl-L-methionine (49), O-phospho-L-serine (37)", "omega-N,omega-N'-dimethyl-L-arginine (67)", ""),
      Tool = "ProSight",
      ConversionList = list("N-acetyl-L-methionine" = "Acetyl")
    ),
    "ConversionList does not have a modification for O-phospho-L-serine"
  )

  expect_message(
    create_proforma(
      Sequence = c("GRGKTGGKARAKAKSRSSRAGLQFPVGRVHRL", "KKTRIIPRHLQLAIRNDEELNKLLGGVTIAY", "TEST"),
      Modifications = c("Methyl 8,Phospho 12,Phospho 22,DiMethyl 30", "TriMethyl 6,Phospho 15", ""),
      Tool = "MSPathFinder",
      ConversionList = list("N-acetyl-L-methionine" = "Acetyl", "O-phospho-L-serine" = "Phospho",
                            "omega-N,omega-N'-dimethyl-L-arginine" = "DiMethyl")
    ),
    "ConversionList is only used with ProSight datasets."
  )
  
  
  ### Running Function ###

  # MSPathFinder example with scan and protein
  MSPathTest <- create_proforma(
    Sequence = c("GRGKTGGKARAKAKSRSSRAGLQFPVGRVHRL", "KKTRIIPRHLQLAIRNDEELNKLLGGVTIAY", "TEST"),
    Modifications = c("Methyl 0,Phospho 12,Phospho 22,DiMethyl 30", "TriMethyl 6,Phospho 15", ""),
    Tool = "MSPathFinder",
    Scan = c(3334, 3336, 3338),
    Protein = c("Protein33", "Protein45", "Protein47")
  )
  expect_true("Proteoform" %in% colnames(MSPathTest))

  # ProSight example without scan and protein
  ProSightTest <- create_proforma(
    Sequence = c("(49)M(37)SGRGKQG", "SG(67)RGKQGGKARAKAKSRSSRAG", "TEST"),
    Modifications = c("N-acetyl-L-methionine (49), O-phospho-L-serine (37)", "omega-N,omega-N'-dimethyl-L-arginine (67)", ""),
    Tool = "ProSight",
    ConversionList = list("N-acetyl-L-methionine" = "Acetyl", "O-phospho-L-serine" = "Phospho",
                          "omega-N,omega-N'-dimethyl-L-arginine" = "DiMethyl")
  )
  expect_true("Proteoform" %in% colnames(ProSightTest))

  # pTop Example
  pTopTest <- create_proforma(
    Sequence = c("SGRGKGGKGLGKGGAKRHRKVLRDNIQGITKPAIRRL", "TEST", "PEPSRSTPAPKKGSKKAITKAQKKDGKKRKRGRKESYSIYV"),
    Modifications = c("(20)Dimethyl[K];(16)Acetyl[K];(0)Acetyl[AnyN-term];", "", "(20)Dimethyl[K];"),
    Tool = "pTop"
  )
  expect_true("Proteoform" %in% colnames(pTopTest))

})
