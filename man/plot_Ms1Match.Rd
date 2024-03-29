% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_Ms1Match.R
\name{plot_Ms1Match}
\alias{plot_Ms1Match}
\title{Plots the Biomolecule Isotope Profile on top of the experimental sequence}
\usage{
plot_Ms1Match(PeakData, Ms1Match, ID, Window = 2)
}
\arguments{
\item{PeakData}{(peak_data object) A pspecterlib peak_data object or data.table with "M/Z" and "Intensity". Required.}

\item{Ms1Match}{(IsoMatchMS_MatchedPeaks) An IsoMatchMS_MatchedPeaks object from match_biomolecule_to_ms1. Required.}

\item{ID}{(numeric) The ID in the IsoMatchMS_MatchedPeaks object to plot. Required.}

\item{Window}{(numeric) The -/+ m/z value on either side of the matched spectra plot. Default is 5 m/z.}
}
\value{
(ggplot object) A ggplot object
}
\description{
Returns a static plot with identified calculated peaks
plotted over the experimental spectrum.
}
\examples{
\dontrun{

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
plot_Ms1Match(
  PeakData = PeakData,
  Ms1Match = IsoMatch,
  ID = 1
)

}

}
