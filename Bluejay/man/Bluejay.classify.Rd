% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/Bluejay_v0.21.R
\name{Bluejay.classify}
\alias{Bluejay.classify}
\title{Classify OV patients drug response by CA125 profile}
\usage{
Bluejay.classify(inputfile, patient_name, Days_plot = 1000, Day_limit = 200,
  CA125_limit = 500, CA125_bound = 35, dfs = 21)
}
\arguments{
\item{inputfile}{CA125 pd file}

\item{patient_name}{name of the patient}

\item{Days_plot}{days after the first CA125 test}

\item{Day_limit}{time period, this function will classify patient base on the CA125 history of these days + 21 days (3 weeks) after surgery}

\item{CA125_limit}{CA125 upper limit when plotting}

\item{CA125_bound}{CA125 value lower bound, we will consider this patient temporary cured if her CA125 value lower than this number}

\item{dfs}{default 21 days (3 weeks) after surgery}
}
\value{
a drug response label: non-determined, sensitive, resistant
}
\description{
Plot CA125 profile and drug response label
}
\details{
This function classifies drug response by using CA125 history profile
}

