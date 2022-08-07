#'@title predict.penCSC
#'
#'@description Flexible prediction method for the objects of class `penCSC` including the absolute risk prediction.
#'
#'@author Shahin Roshani
#'
#'@param object An object of class `penCSC`.
#'@param newX A data frame containing the information of variables related to new records. Information of variables not included in the model creation will be ignored.
#'@param event A vector of event codes which we want predictions for. This must be the same as values in the status variable of the data that was used to create the models. If \code{NULL}, absolute risk will be calculated for all involved events. Default is \code{NULL} which returns values for all involved causes.
#'@param time A vector of time horizons which we want absolute risk predictions at. Only applicable when \code{type='absRisk'}.
#'@param type Type of the predictions. Valid values are: \code{'lp'} or \code{'link'} for linear predictors, \code{'risk'} or \code{'response'} for \code{exp(lp)} and finally \code{'absRisk'} for semi-parametric estimates of absolute risk.
#'@param reference Reference for centering predictions. Valid values are \code{'zero'} and \code{'sample'}. Default is \code{'zero'}. For more information on referencing see details in \code{?predict.coxph}.
#'@param ... Additional arguments. Not used by \code{predict.penCSC}.
#'
#'@return A tibble containing the predictions based on the input arguments.
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
#'predict(penfit,Melanoma[1:5,],type='lp')
#'
#'predict(penfit,Melanoma[1:5,],type='response')
#'
#'predict(penfit,Melanoma[1:5,],type='absRisk',event=1:2,time=1825*(1:2))
#'
#'}
#'
#'@references Pfeiffer, R. M., & Gail, M. M. (2017). Absolute risk: Methods and applications in clinical management and public health.
#'
#'Aalen, O.O. (1978) Nonparametric Inference for a Family of Counting Processes. The Annals of Statistics, 6, 701-726. \doi{10.1214/aos/1176344247}.
#'
#'Wickham H, Averick M, Bryan J, Chang W, McGowan L, François R, et al. Welcome to the tidyverse. J Open Source Softw. 2019 Nov 21;4(43):1686.
#'
#'Bache S, Wickham H (2022). magrittr: A Forward-Pipe Operator for R. \url{https://magrittr.tidyverse.org}, \url{https://github.com/tidyverse/magrittr}.
#'
#'Friedman J, Hastie T, Tibshirani R (2010). "Regularization Paths for Generalized Linear Models via Coordinate Descent." Journal of Statistical Software, 33(1), 1-22. \doi{10.18637/jss.v033.i01}, \url{https://www.jstatsoft.org/v33/i01/}.
#'
#'@import tidyverse survival riskRegression prodlim magrittr glmnet
#'
#'@export

predict.penCSC <- function(object,newX,event=NULL,time,

                           type='lp',reference='zero',...){

  if (!(type %in% c('lp','risk','link','response','absRisk'))){

    stop('type must be `lp`/`link`, `risk`/`response` or `absRisk`!',call.=F)

  }

  if (is.null(event)) event <- names(object$models) %>% stringr::str_remove('Event: ')

  stopifnot('reference must be either `sample` or `zero`!'=reference %in% c('sample','zero'))

  #if (!all(names(newX) %in% names(object$data$input.data))){

  #  stop('newX must have all predictors used in the model!',call.=F)

  #}

  newX <- tibble::as_tibble(newX) %>% dplyr::mutate_if(is.character,as.factor)

  fcts <- object$data$input.data %>% dplyr::select_if(~!is.numeric(.)) %>% purrr::map(levels)

  if (!purrr::is_empty(fcts)){

    for (i in seq_len(length(fcts))){

      levels(newX[[names(fcts)[i]]]) <- fcts[[names(fcts)[i]]]

    }

  }


  lp_risk_pred <- function(object,newX,event,type,reference){


    newX <- object$predictors[event] %>% purrr::map(~model.matrix(.,data=newX)[,-1])


    if (type=='lp') type <- 'link' ; if (type=='risk') type <- 'response'

    if (is.null(event)) event <- names(object$models) %>% stringr::str_remove('Event: ')

    preds <- purrr::pmap(.l=list(object$models[stringr::str_c('Event: ',event)],

                                 newX,

                                 object$parameters$lambda.list[event]),

                         .f=~predict(..1,newx=..2,s=..3,type='link')) %>%

      purrr::map(~tibble::as_tibble(.) %>% dplyr::rename('prediction'='1'))


    if (reference=='sample'){

      means <- object$data$X[event] %>% purrr::map(function(x){

        x <- as.data.frame(x)

        x %>% purrr::map_if(~all(unique(na.omit(.)) %in% 0:1),~0,.else=mean) %>%

          as.data.frame %>% unlist

      })

      pred_modif <- purrr::map2(.x=object$coefs[stringr::str_c('Event: ',event)] %>%

                                  purrr::map(~.[,1] %>% as.vector),

                                .y=means,

                                .f=~sum(.x*.y,na.rm=TRUE))


      preds <-  purrr::map2(.x=preds,

                            .y=pred_modif,

                            .f=~dplyr::mutate_all(.x,function(a) a-.y))

    }


    if (type=='response') preds <- preds %>% purrr::map(~dplyr::mutate_all(.,exp))

    purrr::map2(.x=preds %>% purrr::map(~dplyr::mutate(.,id=seq_len(nrow(.)))),

                .y=names(preds),.f=~dplyr::mutate(.x,event=.y %>% stringr::str_remove('Event: '))) %>%

      purrr::reduce(rbind) %>% dplyr::relocate('prediction',.after='event') -> preds

    return(preds)

  }

  if (type=='absRisk'){

    CSCabsRisk <- function(object,newdata,events,horizons){

      stopifnot('object must be from class penCSC!'=inherits(object,'penCSC'))

      indivs_absRisk <- function(object,indivs,event,horizon){

        dat <- object$data$input.data

        event_times <- dat[which(dat[[object$surv.names[['time']]]]<=horizon &

                                   dat[[object$surv.names[['event']]]]==event),] %>%

          (function(x) x[[object$surv.names[['time']]]])


        b <- object$baseline_hazards[[stringr::str_c('Event: ',event)]] %>%

          dplyr::filter(time %in% event_times) %>% (function(x) x$haz)


        purrr::map2(.x=object$predictors,

                    .y=rep(list(indivs),length(object$predictors)),

                    .f=~model.matrix(.x,data=.y)[,-1]) -> ndX


        cumulative_hazards <- purrr::pmap(.l=list(object$models,

                                                  object$data$X,

                                                  object$data$y,

                                                  ndX,

                                                  object$parameters$lambda.list),

                                          .f=~survival::survfit(..1,

                                                                x=..2,

                                                                y=survival::Surv(time=..3[,1],event=..3[,2]),

                                                                newx=..4,

                                                                s=..5)) %>%

          purrr::map((function(x){

            tibble::tibble(time=x$time,cumhaz=x$cumhaz) %>% dplyr::filter(time %in% event_times)

          })) %>% purrr::map(~.$cumhaz %>% as.matrix)


        risks <- lp_risk_pred(object,indivs,event,'risk',reference='zero') %>% (function(x) x$prediction)


        risks * (matrix(b,nrow=1) %*% (purrr::reduce(cumulative_hazards,`+`) %>%

                                         (function(x) apply(x,2,function(y) exp(-y))))) %>% as.vector -> absRisks


        return(

          tibble::tibble(id=seq_len(nrow(indivs)),event=event,horizon=horizon,absoluteRisk=absRisks)

        )

      }

      grid <- expand.grid(event=events,horizon=horizons) %>% dplyr::mutate_if(is.factor,as.character)

      return(grid %>% purrr::pmap(~indivs_absRisk(object,newdata,..1,..2)) %>%

               purrr::reduce(rbind) %>% tibble::as_tibble())

    }

    preds_res <- CSCabsRisk(object,newX,event,time)

  } else{

    preds_res <- lp_risk_pred(object,newX,event,type,reference)

  }

  return(preds_res)

}

