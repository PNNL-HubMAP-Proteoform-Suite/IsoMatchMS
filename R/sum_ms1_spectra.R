#' Sum MS1 Spectra
#' 
#' @description Reads an mzML file, extracts all MS1 spectra, and sums the abundances
#'     together at a PPM threshold
#' 
#' @param mzMLPath The path to the mzML file. Required.
#' @param PPM_Round Number to round to in PPM. Default is 5. 
#' @param MinimumAbundance The minimum abundance value to keep. Acceptable values
#'    range from 0 to 100. Default is 0.01. 
#' @param Messages A TRUE/FALSE to indicate whether help messages should be printed or not. 
#'    Default is TRUE.
#'    
#' @details
#' Allows users to easily sum spectra for the main IsoMatchMS pipeline function 
#' run_isomatchms. Values are rounded to the nearest 5 ppm. After summing, abundances
#' less than 0.01 are tossed.
#' 
#' @returns A peak_data object with summed MS1 spectra. This object can be directly
#' inputted into run_isomatchms.
#' 
#' @examples
#' \dontrun{
#' 
#' 
#' }
#' @export
sum_ms1_spectra <- function(mzMLPath,
                            PPM_Round = 5,
                            MinimumAbundance = 0.01,
                            Messages = FALSE) {
  
  ##################
  ## CHECK INPUTS ##
  ##################
  
  # pspecterlib will check the mzML file
  
  # PPM Round should not be a value that will slow down the functions too much
  if (!is.numeric(PPM_Round) | PPM_Round < 0.1) {
    stop("PPM_Round must be a numeric of length 1")
  }
  
  # Check that the minimum abundance filter is reasonable 
  if (!is.numeric(MinimumAbundance) || (MinimumAbundance < 0 | MinimumAbundance > 100)) {
    stop("MinimumAbundance should be a single numeric between 0 and 100.")
  }
  
  # Assert that Messages is a True or False statement
  if (!is.logical(Messages) || (Messages != TRUE & Messages != FALSE)) {
    stop("Messages must be TRUE or FALSE")
  }
  
  #####################
  ## SUM MS1 Spectra ##
  #####################
  
  # Pull header information
  header <- pspecterlib::get_scan_metadata(MSPath = mzMLPath)
  
  # Pull MS Level 1
  MS1_Scans <- header[header$`MS Level` == 1, "Scan Number"] %>% unlist() %>% sort()
  
  # If there are no MS1_Scans, return NULL
  if (length(MS1_Scans) == 0) {
    warning("No MS1 scan numbers were found.")
    return(NULL)
  }
  
  if (Messages) {message(paste(length(MS1_Scans), "MS1 scans were detected and are being summed."))}
  
  # Function to sum unction to sum to the ppm value
  ppm_round <- function(vals) {
    
    lapply(vals, function(x) {
      
      if (x >= 1e7) {stop(paste0("Unexpectedly large number passed to rounding function: ", x))}
      
      # Get the number of decimal points
      deci_num <- (x %>% as.integer() %>% nchar()) - 1
      
      # Now use plyr's round any
      adj <- plyr::round_any(x, PPM_Round / 10^6 * 10^deci_num)
      
      return(adj)
      
    }) %>% unlist()
    
  }
  
  if (Messages) {message("Pulling peaks...")}
  
  # Iterate through MS1_Scans, pull peak data, and round MZ values
  All_MS1_Rounded <- do.call(rbind, lapply(MS1_Scans, function(MS1) {
    if (Messages) {
      pos <- match(MS1, MS1_Scans) 
      if (pos %% 500 == 0) {message(paste("...finished", pos, "scans"))}
    }
    charge <- unlist(header[header$`Scan Number` == MS1, "Precursor Charge"])
    peakData <- pspecterlib::get_peak_data(header, MS1, MinAbundance = 0)
    class(peakData) <- c("data.table", "data.frame")
    browser()
    peakData <- peakData %>% 
      dplyr::mutate(
        `M/Z` = (`M/Z` + (charge * 1.00727647)) * charge,
        `M/Z` = ppm_round(`M/Z`)
      ) %>%
      dplyr::select(`M/Z`, Abundance)
    return(peakData)
  }))
  
  if (Messages) {message("Summing peaks...")}
  
  # Sum MS1 scans based on the reounded values
  Sum_MS1 <- All_MS1_Rounded %>%
    dplyr::group_by(`M/Z`) %>%
    dplyr::summarise(
      SumAbundance = sum(Abundance)
    ) %>%
    dplyr::mutate(
      ScaleAbundance = SumAbundance / max(SumAbundance) * 100
    ) %>%
    dplyr::filter(ScaleAbundance > MinimumAbundance)
  
  if (Messages) {message("Peak summing complete!")}
  
  # Generate a peak data object
  return(
    pspecterlib::make_peak_data(
      MZ = Sum_MS1$`M/Z`,
      Intensity = Sum_MS1$ScaleAbundance
    )
  )
  
}