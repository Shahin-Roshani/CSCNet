#'@title predictRisk.penCSC
#'
#'@description predictRisk method for absolute risk prediction. This is mainly for compatibility of 'CSCNet' with functions of 'riskRegression' package.
#'
#'@author Shahin Roshani
#'
#'@param object An object of class 'penCSC'.
#'@param newdata A data frame containing the variable information of new records.
#'@param times A vector of time horizons which we want the absolute risk predictions at.
#'@param cause A single value indicating the event of interest which we want the absolute risk predictions for. This value should be one of the values in the status variable of the data.
#'@param ... Additional arguments. Not used by \code{predictRisk.penCSC}.
#'
#'@return A matrix with columns of absolute risk predictions of individuals for each requested time horizon.
#'
#'@examples \dontrun{
#'
#'data(Melanoma)
#'
#'vl <- list('1'=c('age','sex','ulcer','thick'),
#'
#'           '2'=~age+sex+epicel+thick+ici)
#'
#'al <- list('1'=0,'2'=.5)
#'
#'ll <- list('1'=.01,'2'=.04)
#'
#'penfit <- penCSC(time='time',status='status',vars.list=vl,
#'
#'                 data=Melanoma,alpha.list=al,lambda.list=ll)
#'
#'predictRisk(penfit,Melanoma[1:5,],times=1825*(1:2),cause=1)
#'
#'}
#'
#'@references Wickham H, Averick M, Bryan J, Chang W, McGowan L, François R, et al. Welcome to the tidyverse. J Open Source Softw. 2019 Nov 21;4(43):1686.
#'
#'Bache S, Wickham H (2022). magrittr: A Forward-Pipe Operator for R. \url{https://magrittr.tidyverse.org}, \url{https://github.com/tidyverse/magrittr}.
#'
#'@seealso \url{https://www.rdocumentation.org/packages/riskRegression/versions/1.3.7/topics/predictRisk}
#'
#'Details in: \url{https://rdrr.io/cran/riskRegression/man/Score.html}
#'
#'@import tidyverse survival riskRegression prodlim magrittr glmnet
#'
#'@export

predictRisk.penCSC <- function(object,newdata,times,cause,...){

  if (length(cause)>1){

    stop('`predictRisk` method only handles one cause at a time on `penCSC` objects. To get the results of more times at one time, use `predict.penCSC` method!',call.=F)

  }

  as.list(times) %>%

    purrr::map(~predict.penCSC(object=object,newX=newdata,event=cause,time=.,

                               type='absRisk')$absoluteRisk %>% as.matrix) %>%

    purrr::reduce(cbind) -> absRisks

  colnames(absRisks) <- times

  return(absRisks)

}
