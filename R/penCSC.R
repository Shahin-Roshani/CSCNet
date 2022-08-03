#'@title penCSC
#'
#'@description Function to fit penalized cause-specific-cox with elastic-net penalty.
#'
#'@author Shahin Roshani
#'
#'@param time A character showing the name of the time variable in the data.
#'@param status A character showing the name of the status/event variable in the data.
#'@param vars.list A named list containing the variables to be included in each cause-specific model. Variables can be vectors of variable names or a one sided formula. Names of the list must be the events and exactly the same as values in the status variable. See `Examples` for details.
#'@param data A data frame containing the information of the variables.
#'@param alpha.list A named list containing the single alpha values of each cause-specific model. Names of the list must be the events and exactly the same as values in the status variable. See `Examples` for details.
#'@param lambda.list A named list containing the single lambda values of each cause-specific model. Names of the list must be the events and exactly the same as values in the status variable. See `Examples` for details.
#'@param standardize Logical indicating whether the variables must be standardized or not. Default is \code{TRUE}.
#'@param keep A character vector of the names of variables that should not be shrunk. Default is \code{NULL}.
#'
#'@return A named list containing all the information related to the used data and the fitted models for all causes. Use \code{$} to explore all the involved information.
#'
#'@examples data(Melanoma)
#'
#'vl <- list('1'=c('age','sex','ulcer','thick'),'2'=~age+sex+epicel+thick+ici)
#'
#'al <- list('1'=0,'2'=.5)
#'
#'ll <- list('1'=.01,'2'=.04)
#'
#'penCSC(time='time',status='status',vars.list=vl,data=Melanoma,alpha.list=al,lambda.list=ll)
#'
#'@references Friedman J, Hastie T, Tibshirani R (2010). "Regularization Paths for Generalized Linear Models via Coordinate Descent." Journal of Statistical Software, 33(1), 1-22. doi: 10.18637/jss.v033.i01, \url{https://www.jstatsoft.org/v33/i01/}.
#'
#'@import tidyverse survival riskRegression prodlim magrittr glmnet
#'
#'@export

penCSC <- function(time,status,vars.list,data,

                   alpha.list,lambda.list,standardize=T,keep=NULL){

  name_error <- function(x,arg){

    if (is.null(names(x)) | any(is.na(names(x))) | any(names(x)=='')){

      stop(str_c('Please specify complete names for ',arg,

                 ' to know which event they are related to!'),call.=F)

    }

  }

  name_error(vars.list,'vars.list')

  name_error(alpha.list,'alpha.list')

  name_error(lambda.list,'lambda.list')


  data <- as_tibble(data) %>% mutate_if(is.character,as.factor)


  cens.ind <- data[[status]] %>% unique %>% .[-which(. %in% names(vars.list))]

  y_mats <- names(vars.list) %>%

    as.list %>% (function(x){names(x) <- x ; return(x)}) %>%

    map(~Surv(data[[time]],data[[status]]==.) %>% as.matrix)


  vars.list <- vars.list %>% map(function(x){

    if (inherits(x,'character')) x <- str_c('~',str_c(x,collapse='+')) %>% as.formula

    return(x)

  })

  X_mats <- vars.list %>% map(~model.matrix(.,data=data)[,-1,drop=F])


  alpha.list <- alpha.list[names(y_mats)]


  lambda.list <- lambda.list[names(y_mats)]


  if (!is.null(keep)){

    keep <- keep[names(y_mats)] %>%

      map_if(function(x) !is_empty(x),~model.matrix(str_c('~',str_c(.,collapse='+')) %>%

                                                      as.formula,data=data) %>% colnames %>% .[-1])

    penalty_factors <- map2(.x=keep,.y=X_mats %>% map(colnames),.f=~as.numeric(!(.y %in% .x)))

  } else{

    penalty_factors <- X_mats %>% map(~rep(1,ncol(.)))

  }


  model_lambdas <- pmap(.l=list(aa=y_mats,bb=X_mats,cc=alpha.list,dd=penalty_factors,ee=lambda.list),

                        .f=function(aa,bb,cc,dd,ee){

                          glmnet(x=bb,y=aa,family='cox',alpha=cc,penalty.factor=dd,

                                 standardize=standardize)$lambda %>% c(.,ee) %>%

                            sort(decreasing=T) %>% unique()

                        })


  fits <- pmap(.l=list(y_mats,X_mats,alpha.list,model_lambdas,penalty_factors),

               .f=~glmnet(x=..2,y=..1,family='cox',alpha=..3,lambda=..4,penalty.factor=..5,

                          standardize=standardize))


  names(fits) <- str_c('Event: ',names(fits))


  baseline_ndX <- vars.list %>% map(~model.matrix(.,data=data[1,])[,-1]*0)

  cumulative_baseline_hazards <-

    pmap(.l=list(fits,X_mats,y_mats,baseline_ndX,lambda.list),

         .f=~survfit(..1,x=..2,y=Surv(time=..3[,1],event=..3[,2]),newx=..4,s=..5)) %>%

    map((function(x){

      tibble(time=x$time,cumhaz=x$cumhaz)

    }))


  baseline_hazards <- cumulative_baseline_hazards %>%

    map(function(x){

      x$cumhaz <- diff(c(0,x$cumhaz))

      names(x)[2] <- 'haz'

      return(x)

    })


  coefs <- pmap(.l=list(fits,lambda.list),.f=~predict(..1,s=..2,type='coefficients'))


  result <- list('data'=list('input.data'=data,'y'=y_mats,'X'=X_mats),

                 'models'=fits,'coefs'=coefs,'predictors'=vars.list,

                 'surv.names'=c('time'=time,'event'=status),

                 'parameters'=list('alpha.list'=alpha.list,'lambda.list'=lambda.list),

                 'baseline_hazards'=baseline_hazards)

  class(result) <- 'penCSC'

  print.penCSC <<- function(x) print(x$coefs)

  return(result)

}
