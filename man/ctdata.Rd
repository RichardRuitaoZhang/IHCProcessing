% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IHCProcessing.R
\name{ctdata}
\alias{ctdata}
\title{Extracting and output labeled count data for downward analysis}
\usage{
ctdata(data, output = FALSE, name = NULL)
}
\arguments{
\item{data}{a \code{list} data output from the function \code{labind}}

\item{output}{if you need to output this data for analysis in other software}

\item{name}{the name of your ooutput file}
}
\value{
A data frame with labeled and binded general data.
}
\description{
This function is used to generate the processed labeled count data read from last step,
which can also output the data file for further Prism analyze.
}
\details{
This function is used to extract and output labeled count data for downward analysis, generate
the processed general labeled data read from last step, which can also output the data file for
further Prism analyze. It will return a data frame with labeled and binded general data, and it's
your choice if output it or not.
}
\author{
Ruitao Zhang
}
