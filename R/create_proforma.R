#' Create a ProForma strings from various identification tools
#'
#' @description Derive proforma strings from outputs of MSPathFinder, ProSight, and
#'     pTop. Note that TopPIC outputs ProForma strings so this step is
#'     not required. For mzid files, use pull_modifications_from_mzid.
#'
#' @param Sequence (character) vector containing the "Sequence" column from MSPathFinder,
#'    ProSight, or pTop. Required.
#' @param Modifications (character) vector containing the "Modifications" column
#'    from MSPathFinder, the "PTMs" column from ProSight, or the "PTMs" column from
#'    pTop. Must be the same length as Sequence. If there is no modification, use a blank string. Required.
#' @param Tool (character) indicates which tool was used. Acceptable inputs are "MSPathFinder".
#'    "ProSight", or "pTop". Required.
#' @param ConversionList (list) used to convert ProSight modifications to
#'    ProForma format. See example code and the pspecter modifications glossary with
#'    read.csv(system.file("extdata", "Unimod_v20220602.csv", package = "pspecterlib")). Only
#'    required for ProSight modifications.
#' @param Scan (numeric) vector containing the scan numbers. Must be the same length
#'    as Sequence. Optional.
#' @param Protein (character) vector containing the protein names. Must be the same
#'    length as Sequence. Optional.
#'
#' @returns (data.frame) Returns a data.frame of scan numbers, proforma strings, and protein names.
#'     Note that all unmodified sequences are included as well.
#'
#' @examples
#' \dontrun{
#'
#' # MSPathFinder example with scan and protein
#' create_proforma(
#'    Sequence = c("GRGKTGGKARAKAKSRSSRAGLQFPVGRVHRL", "KKTRIIPRHLQLAIRNDEELNKLLGGVTIAY", "TEST"),
#'    Modifications = c("Methyl 8,Phospho 12,Phospho 22,DiMethyl 30", "TriMethyl 6,Phospho 15", ""),
#'    Tool = "MSPathFinder",
#'    Scan = c(3334, 3336, 3338),
#'    Protein = c("Protein33", "Protein45", "Protein47")
#' )
#'
#' # ProSight example without scan and protein
#' create_proforma(
#'    Sequence = c("(49)M(37)SGRGKQG", "SG(67)RGKQGGKARAKAKSRSSRAG", "TEST"),
#'    Modifications = c("N-acetyl-L-methionine (49), O-phospho-L-serine (37)", "omega-N,omega-N'-dimethyl-L-arginine (67)", ""),
#'    Tool = "ProSight",
#'    ConversionList = list("N-acetyl-L-methionine" = "Acetyl", "O-phospho-L-serine" = "Phospho",
#'                          "omega-N,omega-N'-dimethyl-L-arginine" = "DiMethyl")
#' )
#'
#' # pTop Example
#' create_proforma(
#'    Sequence = c("SGRGKGGKGLGKGGAKRHRKVLRDNIQGITKPAIRRL", "TEST", "PEPSRSTPAPKKGSKKAITKAQKKDGKKRKRGRKESYSIYV"),
#'    Modifications = c("(20)Dimethyl[K];(16)Acetyl[K];(0)Acetyl[AnyN-term];", "", "(20)Dimethyl[K];"),
#'    Tool = "pTop"
#' )
#'
#' }
#' @export
create_proforma <- function(Sequence,
                            Modifications,
                            Tool,
                            ConversionList = NULL,
                            Scan = NULL,
                            Protein = NULL) {

  ##################
  ## CHECK INPUTS ##
  ##################

  # Sequence should be a vector of characters
  if (!is.character(Sequence)) {
    stop("Sequence should be a vector of character strings.")
  }

  # Sequence should all be valid sequences, with the exception of ProSight.
  SequenceValidity <- lapply(Sequence, pspecterlib::is_sequence) %>% unlist()
  if (!all(SequenceValidity) & Tool != "ProSight") {
    stop(paste0("Not all sequences in Sequence are valid. The first example is: ",
       Sequence[which(SequenceValidity == FALSE)[1]], ". See ?pspecterlib::is_sequence for more details."
    ))
  }

  # Modifications should be a vector of characters
  if (!is.character(Modifications)) {
    stop("Modifications must be a vector of character strings.")
  }

  # Modifications and Sequence should be the same length
  if (length(Sequence) != length(Modifications)) {
    stop("Sequence and Modifications should be the same length.")
  }

  # Tool must be "MSPathFinder", "ProSight", or "pTop"
  if (!(Tool %in% c("MSPathFinder", "ProSight", "pTop"))) {
    stop(paste(Tool, "is not a recognized type of Tool. If you would like to add it, please make an issue on the github page."))
  }

  # ConversionList is only required if the tool is ProSight
  if (Tool == "ProSight") {
    if (is.null(ConversionList)) {
      stop("ConversionList is required for ProSight datasets.")
    }
  } else {
    if (!is.null(ConversionList)) {
      message("ConversionList is only used with ProSight datasets.")
    }
  }

  # Scan is not a requirement
  if (!is.null(Scan)) {

    # Scan numbers should be a be numeric
    if (!is.numeric(Scan)) {
      stop("Scan should be a vector of numeric scan numbers.")
    }

    # Scan should be the same length as sequence
    if (length(Sequence) != length(Scan)) {
      stop("Sequence and Scan should be the same length.")
    }

  }

  # Proteins are not a requirement
  if (!is.null(Protein)) {

    # Scan numbers should be a be numeric
    if (!is.character(Protein)) {
      stop("Protein should be a vector of characters.")
    }

    # Scan should be the same length as sequence
    if (length(Sequence) != length(Protein)) {
      stop("Sequence and Protein should be the same length.")
    }

  }

  #############################
  ## CREATE PROFORMA STRINGS ##
  #############################

  # Set Scan and Protein vectors to NA if they're not shared
  if (is.null(Scan)) {Scan <- rep(NA, length(Sequence))}
  if (is.null(Protein)) {Protein <- rep(NA, length(Sequence))}

  if (Tool == "MSPathFinder") {

    # Apply the modifications to the sequences
    PFs <- lapply(1:length(Sequence), function(pos) {

      # If the modification is blank, return the sequence as is
      if (Modifications[pos] == "") {
        return(Sequence[pos])
      } else {

        # Pull out all the modifications
        theMods <- Modifications[pos] %>% strsplit(",") %>% unlist()

        # Split the sequence
        SeqSplit <- Sequence[pos] %>% strsplit("") %>% unlist()

        for (mod in theMods) {

          # Split the modification from the numeric
          ModSplit <- mod %>% strsplit(" ") %>% unlist()

          # If the residue is 0, paste before position 1. Otherwise, fill after appropriate position
          if (as.numeric(ModSplit[2]) == 0) {
            SeqSplit[1] <- paste0("[", ModSplit[1], "]", SeqSplit[1])
          } else {
            SeqSplit[as.numeric(ModSplit[2])] <- paste0(SeqSplit[as.numeric(ModSplit[2])], "[", ModSplit[1], "]")
          }

        }

        # Return ProForma string
        return(paste0(SeqSplit, collapse = ""))

      }

    }) %>% unlist()

  } else if (Tool == "ProSight") {

    PFs <- lapply(1:length(Sequence), function(pos) {

      # If the modification is blank, return the sequence as is
      if (Modifications[pos] == "") {
        return(Sequence[pos])
      } else {

        # Split the modifications
        ModSplit <- Modifications[pos] %>% strsplit(", ") %>% unlist()

        # Apply modifications
        for (mod in ModSplit) {

          # Get the key
          key <- mod %>% strsplit(" ") %>% unlist() %>% tail(1)

          # Convert to the more common modification name
          mod_name <- ConversionList[strsplit(mod, " ") %>% unlist() %>% head(1)]

          # Check that the conversion list has this compound
          if (mod_name %>% unlist() %>% is.null()) {
            stop(paste0("ConversionList does not have a modification for ", strsplit(mod, " ") %>% unlist() %>% head(1) ))
          }

          # Apply the modification
          Sequence[pos] <- gsub(key, paste0("[", mod_name, "]"), Sequence[pos], fixed = T)

        }

        return(Sequence[pos])

      }

    }) %>% unlist()

  } else if (Tool == "pTop") {

    PFs <- lapply(1:length(Sequence), function(pos) {

      # If the modification is blank, return the sequence as is
      if (Modifications[pos] == "") {
        return(Sequence[pos])
      } else {

        # Remove bracketed terms and split the modifications
        ModSplit <- gsub("\\[.*?\\]", "", Modifications[pos]) %>% strsplit(";") %>% unlist()

        # Split Sequence
        SeqSplit <- strsplit(Sequence[pos], "") %>% unlist()

        # Get modification positions
        for (x in ModSplit) {

          # Pull the residue and the modification name
          residue <- gsub("[^[:digit:]]", "", x) %>% as.numeric()
          mod <- gsub("[^[:alpha:]]", "", x)

          # If the residue is 0, paste before position 1
          if (residue == 0) {
            SeqSplit[1] <- paste0("[", mod, "]", SeqSplit[1])
          } else {
            SeqSplit[residue] <- paste0(SeqSplit[residue], "[", mod, "]")
          }

        }

        # Return ProForma string
        return(paste0(SeqSplit, collapse = ""))

      }
    }) %>% unlist()

  }

  # Return ProForma dataframe
  return(
      data.frame(
        Scan = Scan,
        Proteoform = PFs,
        Protein = Protein
      )
  )

}
