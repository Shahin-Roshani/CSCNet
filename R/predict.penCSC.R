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
#'@examples data(Melanoma)
#'
#'vl <- list('1'=c('age','sex','ulcer','thick'),'2'=~age+sex+epicel+thick+ici)
#'
#'al <- list('1'=0,'2'=.5)
#'
#'ll <- list('1'=.01,'2'=.04)
#'
#'penfit <- penCSC(time='time',status='status',vars.list=vl,data=Melanoma,alpha.list=al,lambda.list=ll)
#'
#'predict(penfit,Melanoma[1:5,],type='lp')
#'
#'predict(penfit,Melanoma[1:5,],type='response')
#'
#'predict(penfit,Melanoma[1:5,],type='absRisk',event=1:2,time=1825*(1:2))
#'
#'@references Pfeiffer, R. M., & Gail, M. M. (2017). Absolute risk: Methods and applications in clinical management and public health.
#'
#'Aalen, O.O. (1978) Nonparametric Inference for a Family of Counting Processes. The Annals of Statistics, 6, 701-726. \url{http://dx.doi.org/10.1214/aos/1176344247}.
#'
#'@import tidyverse survival riskRegression prodlim magrittr glmnet
#'
#'@export

predict.penCSC <- function(object,newX,event=NULL,time,

                           type='lp',reference='zero',...){

  if (!(type %in% c('lp','risk','link','response','absRisk'))){

    stop('type must be `lp`/`link`, `risk`/`response` or `absRisk`!',call.=F)

  }

  if (is.null(event)) event <- names(object$models) %>% str_remove('Event: ')

  stopifnot('reference must be either `sample` or `zero`!'=reference %in% c('sample','zero'))

  #if (!all(names(newX) %in% names(object$data$input.data))){

  #  stop('newX must have all predictors used in the model!',call.=F)

  #}

  newX <- as_tibble(newX) %>% mutate_if(is.character,as.factor)

  fcts <- object$data$input.data %>% select(where(~!is.numeric(.))) %>% map(levels)

  if (!is_empty(fcts)){

    for (i in seq_len(length(fcts))){

      levels(newX[[names(fcts)[i]]]) <- fcts[[names(fcts)[i]]]

    }

  }


  lp_risk_pred <- function(object,newX,event,type,reference){


    newX <- object$predictors[event] %>% map(~model.matrix(.,data=newX)[,-1])


    if (type=='lp') type <- 'link' ; if (type=='risk') type <- 'response'

    if (is.null(event)) event <- names(object$models) %>% str_remove('Event: ')

    preds <- pmap(.l=list(object$models[str_c('Event: ',event)],

                          newX,

                          object$parameters$lambda.list[event]),

                  .f=~predict(..1,newx=..2,s=..3,type='link')) %>%

      map(~as_tibble(.) %>% rename('prediction'='1'))


    if (reference=='sample'){

      means <- object$data$X[event] %>% map(function(x){

        x <- as.data.frame(x)

        x %>% map_if(~all(unique(na.omit(.)) %in% 0:1),~0,.else=mean) %>%

          as.data.frame %>% unlist

      })

      pred_modif <- map2(.x=object$coefs[str_c('Event: ',event)] %>% map(~.[,1] %>% as.vector),

                         .y=means,

                         .f=~sum(.x*.y,na.rm=TRUE))


      preds <-  map2(.x=preds,

                     .y=pred_modif,

                     .f=~mutate_all(.x,function(a) a-.y))

    }


    if (type=='response') preds <- preds %>% map(~mutate_all(.,exp))

    map2(.x=preds %>% map(~mutate(.,id=seq_len(nrow(.)))),

         .y=names(preds),.f=~mutate(.x,event=.y %>% str_remove('Event: '))) %>%

      reduce(rbind) %>% relocate('prediction',.after='event') -> preds

    return(preds)

  }

  if (type=='absRisk'){

    CSCabsRisk <- function(object,newdata,events,horizons){

      stopifnot('object must be from class penCSC!'=inherits(object,'penCSC'))

      indivs_absRisk <- function(object,indivs,event,horizon){

        dat <- object$data$input.data

        event_times <- dat[which(dat[[object$surv.names[['time']]]]<=horizon &

                                   dat[[object$surv.names[['event']]]]==event),] %>%

          .[[object$surv.names[['time']]]]


        b <- object$baseline_hazards[[str_c('Event: ',event)]] %>%

          filter(time %in% event_times) %>% .$haz


        map2(.x=object$predictors,

             .y=rep(list(indivs),length(object$predictors)),

             .f=~model.matrix(.x,data=.y)[,-1]) -> ndX


        cumulative_hazards <- pmap(.l=list(object$models,

                                           object$data$X,

                                           object$data$y,

                                           ndX,

                                           object$parameters$lambda.list),

                                   .f=~survfit(..1,

                                               x=..2,

                                               y=Surv(time=..3[,1],event=..3[,2]),

                                               newx=..4,

                                               s=..5)) %>%

          map((function(x){

            tibble(time=x$time,cumhaz=x$cumhaz) %>% filter(time %in% event_times)

          })) %>% map(~.$cumhaz %>% as.matrix)


        risks <- lp_risk_pred(object,indivs,event,'risk',reference='zero') %>%

          .$prediction


        risks * (matrix(b,nrow=1) %*% (reduce(cumulative_hazards,`+`) %>%

                                         apply(.,2,function(x) exp(-x)))) %>%

          as.vector -> absRisks


        return(

          tibble(id=seq_len(nrow(indivs)),event=event,horizon=horizon,absoluteRisk=absRisks)

        )

      }

      grid <- expand.grid(event=events,horizon=horizons) %>% mutate_if(is.factor,as.character)

      return(grid %>% pmap(~indivs_absRisk(object,newdata,..1,..2)) %>%

               reduce(rbind) %>% as_tibble)

    }

    preds_res <- CSCabsRisk(object,newX,event,time)

  } else{

    preds_res <- lp_risk_pred(object,newX,event,type,reference)

  }

  return(preds_res)

}

