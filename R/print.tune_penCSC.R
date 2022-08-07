#'@title print.tune_penCSC
#'
#'@description Internal method for printing the objects of class \code{tune_penCSC}.
#'
#'@author Shahin Roshani
#'
#'@param x An object of class \code{tune_penCSC}.
#'@param ... Other arguments. Not used by \code{print.tune_penCSC}.
#'
#'@return A modified print of \code{tune_penCSC} objects.
#'
#'@export

print.tune_penCSC <- function(x,...){

  print(x$final_fits)

}




