#'@title tune_penCSC
#'
#'@description A flexible function for tuning the involved hyper-parameters of a penalized cause-specific-cox model with elastic net penalty using the linking idea.
#'
#'@author Shahin Roshani
#'
#'@param time A character showing the name of the time variable in the data.
#'@param status A character showing the name of the status/event variable in the data.
#'@param vars.list A named list containing the variables to be included in each cause-specific model. Variables can be vectors of variable names or a one sided formula. Names of the list must be the events and exactly the same as values in the status variable. See `Examples` for details.
#'@param data A data frame containing the information of the variables.
#'@param horizons A vector of time horizons which we want the absolute risk predictions to be evaluated at.
#'@param event The value for event of interest which we want the absolute risk predictions to be evaluated for. This must be one of the values in the status variable of the data.
#'@param rhs A right hand sided formula indicating the variables to be used in estimating the inverse probability of censoring weighting (IPCW) model. Default is \code{~1}.
#'@param method Resampling method to be used for hyper-parameter tuning. Values can be: \code{'cv'} for cross validation, \code{'repcv'} for repeated cross validation, \code{'lgocv'} for monte-carlo cross validation, \code{'loocv'} for leave one out cross validation and \code{'boot'} for bootstrap. Default is \code{'cv'}.
#'@param k Number of folds. Only applicable for \code{method='cv'} and \code{method='repcv'}. Default is 10.
#'@param times Repeat number of the resampling process. Only applicable for \code{method='repcv'}, \code{method='lgocv'} and \code{method='boot'}. Default is 25.
#'@param p The fraction of data to be used as the training set during resampling. Only applicable for \code{method='lgocv'}. Default is 0.7.
#'@param metrics Evaluation metric (loss function) to be used. Values can be \code{'Brier'} for IPCW brier score, \code{'AUC'} for IPCW AUC or a vector of both. Default is \code{'Brier'}.
#'@param final.metric The evaluation metric to decide the best hyper-parameters set for the final fits on the whole data. When \code{NULL} which is the default value, it takes the value from \code{metrics}. If both \code{'Brier'} and \code{'AUC'} were specified in metrics and \code{final.metric} is \code{NULL}, \code{'Brier'} will be used.
#'@param alpha.grid A named list containing a sequence of alpha values to be evaluated for each cause-specific model. Names of the list must be the events and exactly the same as values in the status variable. Default is \code{NULL} which orders the function to set \code{seq(0,1,.5)} for all cause-specific models. See `Details` for more information.
#'@param lambda.grid A named list containing a sequence of lambda values to be evaluated for each cause-specific model. Names of the list must be the events and exactly the same as values in the status variable. Default is \code{NULL} which orders the function to calculate exclusive lambda sequences for all causes. See `Details` for more information.
#'@param nlambdas.list A names list of single integers indicating the length of lambda sequences which are calculated automatically by the function for each cause. Only applicable when \code{lambda.grid=NULL}. Default is NULL which sets all lengths to 5. See `Details` for more information.
#'@param grow.by Difference between the values in the growing sequence of lambda values to find the maximum value that makes the null model. Only applicable when \code{lambda.grid=NULL}. Default is 0.01. See `Details` for more information.
#'@param standardize Logical indicating whether the variables must be standardized or not during model fitting procedures. Default is \code{TRUE}.
#'@param keep A character vector of the names of variables that should not be shrunk in all model fitting procedures. Default is \code{NULL}.
#'@param preProc.fun A function that accepts a data frame and returns a modified data that has gone through the user's desired pre-processing steps. All modifications from this function will be done during the resampling procedures to avoid data leakage. Default is \code{function(x) x}.
#'@param parallel Logical indicating whether the tuning process should be performed in parallel or not. Default is \code{FALSE}.
#'@param preProc.pkgs A character vector containing the names of packages that was used in creating user's \code{preProc.fun} while using parallel computation. Only applicable if \code{parallel=T} and \code{preProc.fun} is a user specified function using functions from other packages. See 'Examples' for details.
#'@param preProc.globals A character vector containing names of objects included in \code{preProc.fun} to be considered as global objects while using parallel computation. The most frequent ones are the names of the user specified pre processing function or functions within this function. Only applicable if \code{parallel=T} and \code{preProc.fun} is a user specified function. See 'Examples' for details.
#'@param core.nums Number of CPU cores to be used for parallel computation. Only applicable if \code{parallel=T}. Default is \code{future::availableCores()/2}.
#'
#'@return A list containing the detailed information of the hyper-parameter tuning and the validation process, best combination of hyper-parameters and the final fits based on the whole data using the best obtained hyper-parameters. Use \code{$} to explore all the involved information.
#'
#'@details \code{tune_penCSC} has the ability to automatically determine the candidate sequences of alpha & lambda values. Setting any of \code{alpha.grid} & \code{lambda.grid} to \code{NULL} will order the function to calculate them automatically.
#'The process of determining the lambda values automatically is by:
#'\enumerate{
#'\item Starting from lambda=0, the algorithm fits LASSO models until finding a lambda value that creates a NULL model where all variables were shrunk to be exactly zero.
#'\item The obtained lambda value will be used as the maximum value of a sequence starting from 0. The length of this sequence is controlled by values in \code{nlambdas.list}.
#'}
#'This will be done for each cause-specific model to create exclusive sequences of lambdas for each of them.
#'
#'@examples data(Melanoma)
#'
#'vl <- list('1'=~age+sex+epicel+ici,'2'=c('age','ulcer','thick','invasion'))
#'
#'al <- list('1'=0,'2'=c(.5,1))
#'
#'#External standardization function for pre-processing with data frame as its input and output
#'
#'library(recipes)
#'
#'std.fun <- function(data){
#'
#'  cont_vars <- data %>% select(where(~is.numeric(.))) %>% names
#'
#'  cont_vars <- cont_vars[-which(cont_vars %in% c('time','status'))]
#'
#'  #External functions from recipes package are being used
#'
#'  recipe(~.,data=data) %>%
#'
#'    step_center(all_of(cont_vars)) %>%
#'
#'    step_scale(all_of(cont_vars)) %>%
#'
#'    prep(training=data) %>% juice
#'
#'}
#'
#'test <- tune_penCSC(time = 'time',status = 'status',vars.list = vl,data = Melanoma,horizons = 1825,
#'
#'                    event = 1,method = 'cv',k = 5,metrics = 'AUC',alpha.grid = al,standardize = FALSE,
#'
#'                    preProc.fun = std.fun,parallel = TRUE,preProc.pkgs = 'recipes')
#'
#'test
#'
#'@references Friedman J, Hastie T, Tibshirani R (2010). "Regularization Paths for Generalized Linear Models via Coordinate Descent." Journal of Statistical Software, 33(1), 1-22. doi: 10.18637/jss.v033.i01, \url{https://www.jstatsoft.org/v33/i01/}.
#'
#'Saadati, M, Beyersmann, J, Kopp-Schneider, A, Benner, A. Prediction accuracy and variable selection for penalized cause-specific hazards models. Biometrical Journal. 2018; 60: 288– 306. \url{https://doi.org/10.1002/bimj.201600242}.
#'
#'Gerds TA, Kattan MW (2021). Medical Risk Prediction Models: With Ties to Machine Learning (1st ed.). Chapman and Hall/CRC. \url{https://doi.org/10.1201/9781138384484.}
#'
#'Pfeiffer, R. M., & Gail, M. M. (2017). Absolute risk: Methods and applications in clinical management and public health.
#'
#'Kuhn, M. (2008). Building Predictive Models in R Using the caret Package. Journal of Statistical Software, 28(5), 1 - 26. doi:\url{http://dx.doi.org/10.18637/jss.v028.i05}.
#'
#'Bengtsson H (2021). “A Unifying Framework for Parallel and Distributed Processing in R using Futures.” The R Journal, 13(2), 208–227. doi:10.32614/RJ-2021-048, \url{https://doi.org/10.32614/RJ-2021-048}.
#'
#'Vaughan D, Dancho M (2022). furrr: Apply Mapping Functions in Parallel using Futures. \url{https://github.com/DavisVaughan/furrr}, \url{https://furrr.futureverse.org/}.
#'
#'@import tidyverse survival riskRegression prodlim magrittr glmnet caret furrr future recipes
#'
#'@export

tune_penCSC <- function(time,status,vars.list,data,horizons,event,rhs=~1,

                        method='cv',k=10,times=25,p=.7,metrics='Brier',

                        final.metric=NULL,alpha.grid=NULL,lambda.grid=NULL,

                        nlambdas.list=NULL,grow.by=.01,standardize=T,keep=NULL,

                        preProc.fun=function(x) x,parallel=F,preProc.pkgs=NULL,

                        preProc.globals=NULL,core.nums=future::availableCores()/2){

  if (!(method %in% c('loocv','lgocv','cv','repcv','boot'))) stop('`method` must be `loocv`, `lgocv`, `cv`, `repcv` or `boot`!',call.=F)

  if (!all(metrics %in% c('Brier','AUC')) | length(metrics)>2) stop('metrics must be either `Brier` or `AUC`. It can also be a vector of both.',call.=F)

  if (is.null(final.metric)){

    if (length(metrics)>1){

      final.metric <- 'Brier'

    } else{

      final.metric <- metrics

    }

  }

  if (!(final.metric %in% c('Brier','AUC')) | length(final.metric)!=1) stop('final.metric should be only one of `Brier` or `AUC`.',call.=F)

  if (length(event)!=1) stop('Only the event of interest must be specified for the tuning process!',call.=F)

  if (is.character(preProc.fun)){

    if (length(preProc.fun)>1){

      stop('Only one name of a unified pre-processing function must be given!',call.=F)

    } else{

      preProc.globals <- c(preProc.globals,preProc.fun)

    }

  }

  resampler <- function(method){

    if (method=='loocv'){

      indices <- seq_len(nrow(data))

      train_index_list <- as.list(indices) %>% map(~indices[-.])

    }

    if (method=='lgocv') train_index_list <- caret::createDataPartition(y=data[[status]],times=times,p=p,list=T)

    if (method=='cv') train_index_list <- caret::createFolds(y=data[[status]],k=k,list=T,returnTrain=T)

    if (method=='repcv') train_index_list <- caret::createMultiFolds(y=data[[status]],k=k,times=times)

    if (method=='boot') train_index_list <- caret::createResample(y=data[[status]],times=times,list=T)


    test_index_list <- train_index_list %>% map(~which(!(seq_len(nrow(data)) %in% .)))

    return(list(train_index_list=train_index_list,test_index_list=test_index_list))

  }

  codes <- unique(data[[status]])

  cens.code <- codes[-which(codes %in% names(vars.list))]

  form <- str_c('Hist(',time,',',status,',cens.code=\'',cens.code,'\'',')',

                as.character.POSIXt(rhs)) %>% as.formula

  if (is_empty(lambda.grid)){

    if (is_empty(nlambdas.list)){

      nlambdas.list <- rep(list(5),length(vars.list))

      names(nlambdas.list) <- names(vars.list)

    }

    nlambdas.list <- nlambdas.list[names(vars.list)]


    dd <- as_tibble(data) %>% mutate_if(is.character,as.factor) %>% na.omit

    ymats <- names(vars.list) %>%

      as.list %>% (function(x){names(x) <- x ; return(x)}) %>%

      map(~Surv(dd[[time]],dd[[status]]==.) %>% as.matrix)

    vl <- vars.list %>% map(function(x){

      if (inherits(x,'character')) x <- str_c('~',str_c(x,collapse='+')) %>% as.formula

      return(x)

    })

    Xmats <- vl %>% map(~model.matrix(.,data=dd)[,-1,drop=F])


    lambda_seq <- function(X,y,nlambdas){

      max_lambda <- 0

      range_fit <- glmnet(x=X,y=y,family='cox',alpha=1,standardize=T) %>%

        predict(.,s=max_lambda,type='coefficients')

      while (any(range_fit[,1]!=0)){

        max_lambda <- max_lambda + grow.by

        range_fit <- glmnet(x=X,y=y,family='cox',alpha=1,standardize=T) %>%

          predict(.,s=max_lambda,type='coefficients')

      }

      return(seq(0,max_lambda,length.out=nlambdas))

    }

    lambda.grid <- pmap(.l=list(ymats,Xmats,nlambdas.list),.f=~lambda_seq(..2,..1,..3) %>% unique)

  }

  if (is_empty(alpha.grid)){

    alpha.grid <- rep(list(seq(0,1,.5)),length(vars.list))

    names(alpha.grid) <- names(vars.list)

  }


  alpha.grid <- alpha.grid[names(vars.list)]

  lambda.grid <- lambda.grid[names(vars.list)]

  names(alpha.grid) <- str_c('alpha_',names(alpha.grid))

  names(lambda.grid) <- str_c('lambda_',names(lambda.grid))


  start <- Sys.time()


  grid <- c(alpha.grid,lambda.grid,horizon=list(horizons)) %>% expand.grid


  nn <- length(lambda.grid)

  zl_indices <- apply(grid,1,function(x) which(x[(1:nn)+nn]==0) %>% as.vector)

  for (i in 1:nrow(grid)){

    if (!is_empty(zl_indices[[i]])) grid[i,zl_indices[[i]]] <- 0

  }

  grid <- distinct(grid)


  calc_grid <- grid %>% mutate(combination=seq_len(nrow(.))) %>%

    split(~.$combination) %>%

    map(~select(.,-combination) %>% (function(x){

      y <- list()

      y$alpha.list <- select(x,starts_with('alpha_')) %>%

        rename_all(~str_remove(.,'alpha_')) %>% as.list

      y$lambda.list <- select(x,starts_with('lambda_')) %>%

        rename_all(~str_remove(.,'lambda_')) %>% as.list

      y$horizon <- select(x,horizon) %>% unlist

      return(y)

    }))


  modeling <- function(alpha_list,lambda_list,horizon){

    resamples <- resampler(method)

    training_list <- resamples$train_index_list %>% map(~data[.,] %>% preProc.fun)

    testing_list <- resamples$test_index_list %>% map(~data[.,] %>% preProc.fun)

    pmap(.l=list(aa=training_list,bb=testing_list),

         .f=possibly(.f=function(aa,bb){

           penCSC(time=time,

                  status=status,

                  vars.list=vars.list,

                  data=aa,

                  alpha.list=alpha_list,

                  lambda.list=lambda_list,

                  standardize=standardize,

                  keep=keep) -> fit

           riskRegression::Score(list(fit),data=bb,formula=form,metrics=metrics,

                                 cause=event,times=horizon,null.model=F) %>% .[metrics] %>%

             map(~.$score %>% .[,3] %>% unlist %>% as.vector) %>% as_tibble

         },

         otherwise=matrix(NA,1,length(metrics)) %>%

           (function(x){colnames(x) <- metrics ; return(as_tibble(x))})

         )

    ) %>% (function(x){

      map2(.x=x,

           .y=as.list(names(x)),

           .f=~mutate(.x,'step'=.y) %>% relocate('step',.before=1))

    }) %>% reduce(rbind)

  }


  if (!parallel){

    calc_grid %>% map(~modeling(.$alpha.list,.$lambda.list,.$horizon)) -> lossfun_vals

  } else{

    pkg_envs <- c('tidyverse','magrittr','survival','riskRegression','prodlim',

                  'Publish','glmnet','caret','furrr') %>%

      c(.,preProc.pkgs) %>% unique

    globals <- c('penCSC','predict.penCSC','vars.list','preProc.fun',

                 'resampler','data','predictRisk.penCSC','keep','modeling') %>%

      c(.,preProc.globals) %>% unique

    future::plan(future::multisession(),workers=core.nums)

    cat('\n')

    calc_grid %>% furrr::future_map(~modeling(.$alpha.list,.$lambda.list,.$horizon),

                                    .options=furrr::furrr_options(packages=pkg_envs,globals=globals,seed=T),

                                    .progress=T) -> lossfun_vals

    cat('\n')

  }


  stop <- Sys.time()

  cat('\nProcess was done in:',(stop-start) %>% as.character.POSIXt %>% str_c(.,'.\n\n'))


  lossfun_vals %>% map(~select(.,-step) %>% summarize_all(mean) %>%

                         rename_all(~str_c('mean.',.))) %>%

    reduce(rbind) %>% cbind(grid,.) -> validation_result


  final_params <- validation_result %>% split(~.$horizon) %>%

    map(function(x){

      if (final.metric=='Brier'){

        res <- filter(x,mean.Brier==min(mean.Brier,na.rm=T))

      } else{

        res <- filter(x,mean.AUC==max(mean.AUC,na.rm=T))

      }

      return(res)

    })


  final_fits <- final_params %>%

    map(function(x){

      al <- x[str_subset(names(x),'alpha_')] %>% as.list

      ll <- x[str_subset(names(x),'lambda_')] %>% as.list

      names(al) <- str_remove(names(al),'alpha_')

      names(ll) <- str_remove(names(ll),'lambda_')

      penCSC(time = time,

             status = status,

             vars.list = vars.list,

             data = data %>% preProc.fun,

             alpha.list = al,

             lambda.list = ll)

    })

  tuning_results <- list(lossfun_vals=lossfun_vals,

                         validation_result=validation_result,

                         final_params=final_params,

                         final_fits=final_fits)

  class(tuning_results) <- 'tune_penCSC'

  print.tune_penCSC <<- function(x) print(x$final_fits)

  return(tuning_results)

}

