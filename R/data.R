#' Boston Dataset
#' 
#' The Boston Housing Dataset, derived from information from the U.S. Census
#' Service.
#' 
#' @format ## `Boston`
#' A data.frame with 506 rows and 14 columns:
#' \describe{
#'  \item{x1-x13}{The covariates in the data}
#'  \item{y}{The response variable}
#' }
#' @source Harrison, D. and Rubinfeld, D.L. `Hedonic prices and the demand for clean air', J. Environ. Economics & Management, vol.5, 81-102, 1978
"Boston"

#' Weather Dataset
#' 
#' A subset of data collected as part of the GES DISC datasets for calibrated
#' brightness temperatures.
#' 
#' @format ## `Weather`
#' A data.frame with 2000 rows and 3 columns:
#' \describe{
#'  \item{theta}{Polar angle}
#'  \item{phi}{Azimuthal angle}
#'  \item{Y}{Response variable}
#' }
#' 
#' @source https://disc.gsfc.nasa.gov/datasets/GPM_1CGPMGMI_07/summary?keywords=10.5067%2FGPM%2FGMI%2FGPM%2F1C%2F07
"Weather"

#' Cylindrical Dataset
#' 
#' Synthetic data generated for cylindrical parameters. The response is given by
#' the function f(z, theta) = sin(z)cos(theta) + z*sin(theta)^2, with added noise.
#' 
#' @format ## `Cylinder`
#' A data.frame with 400 rows and 3 columns:
#' \describe{
#'  \item{z}{Height}
#'  \item{theta}{Polar angle}
#'  \item{Y}{Response}
#' }
#' 
"Cylinder"
