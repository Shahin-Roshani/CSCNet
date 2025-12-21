#'@title print.penCSC
#'
#'@description Internal method for printing the objects of class \code{penCSC}.
#'
#'@author Shahin Roshani
#'
#'@param x An object of class \code{penCSC}.
#'@param ... Other arguments. Not used by \code{print.penCSC}.
#'
#'@return A modified print of \code{penCSC} objects.
#'
#'@export

print.penCSC <- function(x,...){

  print(x$coefs)

}
