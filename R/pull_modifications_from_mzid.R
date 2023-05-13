#' Pull proforma strings
#'
#' @description Pull modifications from MS-GF+ and MSPathFinder files (mzid).
#'
#' @param IDPath (character) The path to the mzid file. Required.
#'
#' @returns (data.frame) Returns a data.frame of scan numbers, proforma strings, and protein names.
#'     Note that all unmodified sequences are included as well.
#'
#' @examples
#' \dontrun{
#'
#' # Create a temporary directory and copy example data there
#' tmpdir <- tempdir()
#'
#' files <- c(
#' "https://raw.githubusercontent.com/EMSL-Computing/PSpecteR/master/pspecter_container/TestFiles/BottomUp/BottomUp.mzid",
#' "https://raw.githubusercontent.com/EMSL-Computing/PSpecteR/master/pspecter_container/TestFiles/TopDown/TopDown.mzid"
#' )
#'
#'for (file in files) {
#'  download.file(file, file.path(tmpdir, tail(unlist(strsplit(file, "/")), 1)))
#' }
#'
#' # Download a MS-GF+ example
#' pull_modifications_from_mzid(file.path(tmpdir, "BottomUp.mzid"))
#'
#' # Download an MSPathFinder example
#' pull_modifications_from_mzid(file.path(tmpdir, "TopDown.mzid"))
#'
#' }
#' @export
pull_modifications_from_mzid <- function(IDPath) {

  ##################
  ## CHECK INPUT ##
  #################

  # Assert that the ID file is an mzid file
  if (grepl(".mzid|.mzID", IDPath) == FALSE) {
    stop("ID file must be an mzid file.")
  }
  
  
  # Assert that the ID path exists
  if (file.exists(IDPath) == FALSE) {
    stop("ID file path must exist.")
  }

  ######################
  ## READ THE ID FILE ##
  ######################
  
  # Open the identification file
  ID <- mzR::openIDfile(IDPath)

  # Pull Peptide-Spectrum-Match data
  PSM <- mzR::psms(ID)

  # Get scan numbers, sequences (proteoforms) and protein
  Proforma_DF <- data.frame(
    Scan = PSM$acquisitionNum,
    Proteoform = PSM$sequence,
    Protein = PSM$DatabaseAccess
  ) %>%
    dplyr::arrange(Scan)

  # Pull modification data.
  MOD <- mzR::modifications(ID)

  if (nrow(MOD) != 0) {

    # Add scan number, convert peptideReferences to proforma strings, and add protein information.
    MOD <- MOD %>%
      dplyr::mutate(
        Scan = gsub("controllerType=0 controllerNumber=1 scan=", "", spectrumID) %>% as.numeric(), # Extract scan numbers
        Proteoform = gsub("Pep_|[0-9]", "", peptideRef), # Remove extra text and numerics from peptide references
        Name = lapply(1:length(Scan), function(pos) {
          theName <- .$name[pos]
          ifelse(is.null(theName) || is.na(theName) || theName != "", theName, .$mass[pos])
        }) %>% unlist()
      ) %>%
      dplyr::select(Scan, Proteoform, Name) %>%
      unique() %>%
      dplyr::mutate(
        Proteoform = lapply(1:length(Scan), function(pos) {
          fix <- gsub("+", paste0("[", .$Name[pos], "]"), .$Proteoform[pos], fixed = T)
          fix <- gsub("-", paste0("[", .$Name[pos], "]"), fix, fixed = T)
          return(fix)
        }) %>% unlist()
      ) %>%
      merge(Proforma_DF[,c("Scan", "Protein")], by = "Scan") %>%
      dplyr::select(-Name) %>%
      unique() %>%
      dplyr::arrange(Scan)

    # Bind results
    Proforma_DF <- rbind(Proforma_DF, MOD) %>%
      dplyr::arrange(Scan)

    # Return results
    return(Proforma_DF)

  } else {
    
    return(Proforma_DF)
    
  }

}
