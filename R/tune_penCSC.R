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
#'@param strat.var A single character indicating name of the strata variable to be used to create resamples. If numerical, groups will be specified based on percentiles. Default is \code{NULL} which considers status variable as a factor and creates the resamples based on different levels of it.
#'@param metrics Evaluation metric (loss function) to be used. Values can be \code{'Brier'} for IPCW brier score, \code{'AUC'} for IPCW AUC or a vector of both. Default is \code{'Brier'}.
#'@param final.metric The evaluation metric to decide the best hyper-parameters set for the final fits on the whole data. When \code{NULL} which is the default value, it takes the value from \code{metrics}. If both \code{'Brier'} and \code{'AUC'} were specified in metrics and \code{final.metric} is \code{NULL}, \code{'Brier'} will be used.
#'@param alpha.grid A named list containing a sequence of alpha values to be evaluated for each cause-specific model. Names of the list must be the events and exactly the same as values in the status variable. Default is \code{NULL} which orders the function to set \code{seq(0,1,.5)} for all cause-specific models. See `Details` for more information.
#'@param lambda.grid A named list containing a sequence of lambda values to be evaluated for each cause-specific model. Names of the list must be the events and exactly the same as values in the status variable. Default is \code{NULL} which orders the function to calculate exclusive lambda sequences for all causes. See `Details` for more information.
#'@param nlambdas.list A names list of single integers indicating the length of lambda sequences which are calculated automatically by the function for each cause. Only applicable when \code{lambda.grid=NULL}. Default is NULL which sets all lengths to 5. See `Details` for more information.
#'@param grow.by Difference between the values in the growing sequence of lambda values to find the maximum value that makes the null model. Only applicable when \code{lambda.grid=NULL}. Default is 0.01. See `Details` for more information.
#'@param standardize Logical indicating whether the variables must be standardized or not during model fitting procedures. Default is \code{TRUE}.
#'@param keep A character vector of the names of variables that should not be shrunk in all model fitting procedures. Default is \code{NULL}.
#'@param preProc.fun A function that accepts a data and returns a modified version of it that has gone through the user's desired pre-processing steps. All modifications from this function will be done during the resampling procedures to avoid data leakage. It will modify all training and test set(s) during the validation unless other argument \code{preProc.fun.test} is specified by user and then it only affects the training set(s). Default is \code{function(x) x}. Also see the description of \code{preProc.fun.test} argument.
#'@param preProc.fun.test A function the exact same characteristics and description as \code{preProc.fun} argument. If user specifies a separate function for \code{preProc.fun.test}, it will only affect test set(s) during validation while the function from \code{preProc.fun} will affect the training set(s). Default is \code{NULL} which means function from \code{preProc.fun} will be used on both training and test set(s) during validation. Also see the description of \code{preProc.fun} argument.
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
#'@examples \donttest{
#'
#'#install.packages('collinear')
#'
#'library(riskRegression)
#'
#'data(Melanoma)
#'
#'vl <- list('1'=~age+sex+epicel+ici,
#'
#'           '2'=c('age','ulcer','thick','invasion'))
#'
#'al <- list('1'=0,'2'=c(.5,1))
#'
#'#External standardization function with data frame as its input and output
#'
#'zvr.fun <- function(data){
#'
#'  zv_vars <- identify_zero_variance_variables(df = data,responses = c('time','status'))
#'
#'  return(data %>% select(-all_of(zv_vars)))
#'
#'}
#'
#'set.seed(233)
#'
#'test <- tune_penCSC(time='time',status='status',vars.list=vl,data=Melanoma,horizons=1095,
#'
#'                    event=1,method='cv',k=3,metrics='AUC',alpha.grid=al,standardize=TRUE,
#'
#'                    preProc.fun=zvr.fun,parallel=TRUE,preProc.pkgs='collinear')
#'
#'test
#'
#'}
#'
#'@references Friedman J, Hastie T, Tibshirani R (2010). "Regularization Paths for Generalized Linear Models via Coordinate Descent." Journal of Statistical Software, 33(1), 1-22. \doi{10.18637/jss.v033.i01}, \url{https://www.jstatsoft.org/v33/i01/}.
#'
#'Saadati, M, Beyersmann, J, Kopp-Schneider, A, Benner, A. Prediction accuracy and variable selection for penalized cause-specific hazards models. Biometrical Journal. 2018; 60: 288– 306. \doi{10.1002/bimj.201600242}.
#'
#'Gerds TA, Kattan MW (2021). Medical Risk Prediction Models: With Ties to Machine Learning (1st ed.). Chapman and Hall/CRC. \doi{10.1201/9781138384484}
#'
#'Pfeiffer, R. M., & Gail, M. M. (2017). Absolute risk: Methods and applications in clinical management and public health.
#'
#'Kuhn, M. (2008). Building Predictive Models in R Using the caret Package. Journal of Statistical Software, 28(5), 1–26. \doi{10.18637/jss.v028.i05}.
#'
#'Bengtsson H (2021). “A Unifying Framework for Parallel and Distributed Processing in R using Futures.” The R Journal, 13(2), 208–227. \doi{10.32614/RJ-2021-048}.
#'
#'Vaughan D, Dancho M (2022). furrr: Apply Mapping Functions in Parallel using Futures. \url{https://github.com/DavisVaughan/furrr}, \url{https://furrr.futureverse.org/}.
#'
#'Therneau T (2022). A Package for Survival Analysis in R. R package version 3.3-1, \url{https://CRAN.R-project.org/package=survival}.
#'
#'Wickham H, Averick M, Bryan J, Chang W, McGowan L, François R, et al. Welcome to the tidyverse. J Open Source Softw. 2019 Nov 21;4(43):1686.
#'
#'Bache S, Wickham H (2022). magrittr: A Forward-Pipe Operator for R. \url{https://magrittr.tidyverse.org}, \url{https://github.com/tidyverse/magrittr}.
#'
#'@import tidyverse survival riskRegression prodlim magrittr glmnet furrr recipes
#'
#'@importFrom caret createDataPartition createFolds createMultiFolds createResample
#'
#'@importFrom future plan availableCores
#'
#'@importFrom stats predict
#'
#'@importFrom prodlim Hist
#'
#'@export

tune_penCSC <- function(time,status,vars.list,data,horizons,event,rhs=~1,

                        method='cv',k=10,times=25,p=.7,strat.var=NULL,

                        metrics='Brier',final.metric=NULL,alpha.grid=NULL,

                        lambda.grid=NULL,nlambdas.list=NULL,grow.by=.01,standardize=TRUE,

                        keep=NULL,preProc.fun=function(x) x,preProc.fun.test=NULL,

                        parallel=FALSE,preProc.pkgs=NULL,preProc.globals=NULL,

                        core.nums=future::availableCores()/2){

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

  if (!is.function(preProc.fun)) stop('`preProc.fun` must be a function!',call.=F)

  if (!purrr::is_empty(preProc.fun.test) & !is.function(preProc.fun.test)){

    stop('`preProc.fun.test` must be a function!',call.=F)

  }

  #if (is.character(preProc.fun)){

  #  if (length(preProc.fun)>1){

  #    stop('Only one name of a unified pre-processing function must be given!',call.=F)

  #  } else{

  #    preProc.globals <- c(preProc.globals,preProc.fun)

  #  }

  #}

  if (purrr::is_empty(preProc.fun.test)) preProc.fun.test <- preProc.fun

  if (purrr::is_empty(strat.var)){

    strat.vec <- data[[status]] %>% as.factor

  } else{

    strat.vec <- data[[strat.var]]

  }

  resampler <- function(method){

    if (method=='loocv'){

      indices <- seq_len(nrow(data))

      train_index_list <- as.list(indices) %>% purrr::map(~indices[-.])

    }

    if (method=='lgocv') train_index_list <- caret::createDataPartition(y=strat.vec,times=times,p=p,list=T)

    if (method=='cv') train_index_list <- caret::createFolds(y=strat.vec,k=k,list=T,returnTrain=T)

    if (method=='repcv') train_index_list <- caret::createMultiFolds(y=strat.vec,k=k,times=times)

    if (method=='boot') train_index_list <- caret::createResample(y=strat.vec,times=times,list=T)


    test_index_list <- train_index_list %>% purrr::map(~which(!(seq_len(nrow(data)) %in% .)))

    return(list(train_index_list=train_index_list,test_index_list=test_index_list))

  }

  codes <- unique(data[[status]])

  cens.code <- codes[-which(codes %in% names(vars.list))]

  if (purrr::is_empty(cens.code)) cens.code <- stats::rnorm(1)

  form <- stringr::str_c('Hist(',time,',',status,',cens.code=\'',cens.code,'\'',')',

                         format(rhs)) %>% stats::as.formula()

  if (purrr::is_empty(lambda.grid)){

    if (purrr::is_empty(nlambdas.list)){

      nlambdas.list <- rep(list(5),length(vars.list))

      names(nlambdas.list) <- names(vars.list)

    }

    nlambdas.list <- nlambdas.list[names(vars.list)]


    dd <- tibble::as_tibble(data) %>% dplyr::mutate_if(is.character,as.factor) %>% stats::na.omit()

    ymats <- names(vars.list) %>%

      as.list %>% (function(x){names(x) <- x ; return(x)}) %>%

      purrr::map(~survival::Surv(dd[[time]],dd[[status]]==.) %>% as.matrix())

    vl <- vars.list %>% purrr::map(function(x){

      if (inherits(x,'character')) x <- stringr::str_c('~',stringr::str_c(x,collapse='+')) %>% stats::as.formula()

      return(x)

    })

    Xmats <- vl %>% purrr::map(~model.matrix(.,data=dd)[,-1,drop=F])


    lambda_seq <- function(X,y,nlambdas){

      max_lambda <- 0

      range_fit <- glmnet::glmnet(x=X,y=y,family='cox',alpha=1,standardize=T) %>%

        (function(x) predict(x,s=max_lambda,type='coefficients'))

      while (any(range_fit[,1]!=0)){

        max_lambda <- max_lambda + grow.by

        range_fit <- glmnet::glmnet(x=X,y=y,family='cox',alpha=1,standardize=T) %>%

          (function(x) predict(x,s=max_lambda,type='coefficients'))

      }

      return(seq(0,max_lambda,length.out=nlambdas))

    }

    lambda.grid <- purrr::pmap(.l=list(ymats,Xmats,nlambdas.list),

                               .f=~lambda_seq(..2,..1,..3) %>% unique)

  }

  if (purrr::is_empty(alpha.grid)){

    alpha.grid <- rep(list(seq(0,1,.5)),length(vars.list))

    names(alpha.grid) <- names(vars.list)

  }


  alpha.grid <- alpha.grid[names(vars.list)]

  lambda.grid <- lambda.grid[names(vars.list)]

  names(alpha.grid) <- stringr::str_c('alpha_',names(alpha.grid))

  names(lambda.grid) <- stringr::str_c('lambda_',names(lambda.grid))


  start <- Sys.time()


  grid <- c(alpha.grid,lambda.grid,horizon=list(horizons)) %>% expand.grid()


  nn <- length(lambda.grid)

  zl_indices <- apply(grid,1,function(x) which(x[(1:nn)+nn]==0) %>% as.vector)

  for (i in seq_len(nrow(grid))){

    if (!purrr::is_empty(zl_indices[[i]])) grid[i,zl_indices[[i]]] <- 0

  }

  grid <- dplyr::distinct(grid)


  calc_grid <- grid %>% dplyr::mutate(combination=seq_len(nrow(grid))) %>%

    (function(x) split(x,x$combination)) %>%

    purrr::map(~select(.,-combination) %>% (function(x){

      y <- list()

      y$alpha.list <- dplyr::select(x,dplyr::starts_with('alpha_')) %>%

        dplyr::rename_all(~stringr::str_remove(.,'alpha_')) %>% as.list

      y$lambda.list <- dplyr::select(x,dplyr::starts_with('lambda_')) %>%

        dplyr::rename_all(~stringr::str_remove(.,'lambda_')) %>% as.list

      y$horizon <- dplyr::select(x,horizon) %>% unlist

      return(y)

    }))


  resamples <- resampler(method)

  training_list <- resamples$train_index_list %>% purrr::map(~data[.,] %>% preProc.fun)

  testing_list <- resamples$test_index_list %>% purrr::map(~data[.,] %>% preProc.fun.test)


  modeling <- function(alpha_list,lambda_list,horizon){

    force(alpha_list)

    force(lambda_list)

    force(horizon)

    purrr::pmap(.l=list(aa=training_list,bb=testing_list),

                .f=purrr::possibly(.f=function(aa,bb){

                  penCSC(time=time,

                         status=status,

                         vars.list=vars.list,

                         data=aa,

                         alpha.list=alpha_list,

                         lambda.list=lambda_list,

                         standardize=standardize,

                         keep=keep) -> fit

                  riskRegression::Score(list(fit),data=bb,formula=form,metrics=metrics,

                                        cause=event,times=horizon,null.model=F) %>% (function(x) x[metrics]) %>%

                    purrr::map(~.$score %>% (function(x) x[,3]) %>% unlist %>% as.vector) %>% tibble::as_tibble()

                },

                otherwise=matrix(NA,1,length(metrics)) %>%

                  (function(x){colnames(x) <- metrics ; return(tibble::as_tibble(x))})

                )

    ) %>% (function(x){

      purrr::map2(.x=x,

                  .y=as.list(names(x)),

                  .f=~dplyr::mutate(.x,'step'=.y) %>% dplyr::relocate('step',.before=1))

    }) %>% purrr::reduce(rbind)

  }


  if (!parallel){

    calc_grid %>% purrr::map(function(x) modeling(x$alpha.list,x$lambda.list,x$horizon)) -> lossfun_vals

  } else{

    pkg_envs <- c('tidyverse','magrittr','survival','riskRegression','prodlim',

                  'Publish','glmnet','caret','furrr') %>%

      (function(x) c(x,preProc.pkgs)) %>% unique()

    globals <- c('penCSC','predict.penCSC','vars.list','preProc.fun','preProc.fun.test',

                 'strat.var','strat.vec','resampler','data','predictRisk.penCSC','keep',

                 'modeling') %>%

      (function(x) c(x,preProc.globals)) %>% unique()

    future::plan(future::multisession(),workers=core.nums)

    calc_grid %>% furrr::future_map(function(x) modeling(x$alpha.list,x$lambda.list,x$horizon),

                                    .options=furrr::furrr_options(packages=pkg_envs,globals=globals,seed=T),

                                    .progress=T) -> lossfun_vals

  }


  stop <- Sys.time()

  message(stringr::str_c('\nProcess was done in ',format(stop-start),'.'))


  lossfun_vals %>% purrr::map(~dplyr::select(.,-step) %>% dplyr::summarize_all(mean) %>%

                                dplyr::rename_all(~stringr::str_c('mean.',.))) %>%

    purrr::reduce(rbind) %>% (function(x) cbind(grid,x)) -> validation_result


  final_params <- validation_result %>% (function(x) split(x,x$horizon)) %>%

    purrr::map(function(x){

      if (final.metric=='Brier'){

        res <- dplyr::filter(x,x$mean.Brier==min(x$mean.Brier,na.rm=T))

      } else{

        res <- dplyr::filter(x,x$mean.AUC==max(x$mean.AUC,na.rm=T))

      }

      return(res)

    })


  final_fits <- final_params %>%

    purrr::map(function(x){

      al <- x[stringr::str_subset(names(x),'alpha_')] %>% as.list()

      ll <- x[stringr::str_subset(names(x),'lambda_')] %>% as.list()

      names(al) <- stringr::str_remove(names(al),'alpha_')

      names(ll) <- stringr::str_remove(names(ll),'lambda_')

      penCSC(time = time,

             status = status,

             vars.list = vars.list,

             data = data %>% preProc.fun,

             alpha.list = al,

             lambda.list = ll,

             keep = keep,

             standardize = standardize)

    })

  tuning_results <- list(lossfun_vals=lossfun_vals,

                         validation_result=validation_result,

                         final_params=final_params,

                         final_fits=final_fits)

  class(tuning_results) <- 'tune_penCSC'

  return(tuning_results)

}

