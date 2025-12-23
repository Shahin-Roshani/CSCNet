# tune_penCSC

A flexible function for tuning the involved hyper-parameters of a
penalized cause-specific-cox model with elastic net penalty using the
linking idea.

## Usage

``` r
tune_penCSC(
  time,
  status,
  vars.list,
  data,
  horizons,
  event,
  rhs = ~1,
  method = "cv",
  k = 10,
  times = 25,
  p = 0.7,
  strat.var = NULL,
  metrics = "Brier",
  final.metric = NULL,
  alpha.grid = NULL,
  lambda.grid = NULL,
  nlambdas.list = NULL,
  grow.by = 0.01,
  standardize = TRUE,
  keep = NULL,
  preProc.fun = function(x) x,
  preProc.fun.test = NULL,
  parallel = FALSE,
  preProc.pkgs = NULL,
  preProc.globals = NULL,
  core.nums = future::availableCores()/2
)
```

## Arguments

- time:

  A character showing the name of the time variable in the data.

- status:

  A character showing the name of the status/event variable in the data.

- vars.list:

  A named list containing the variables to be included in each
  cause-specific model. Variables can be vectors of variable names or a
  one sided formula. Names of the list must be the events and exactly
  the same as values in the status variable. See \`Examples\` for
  details.

- data:

  A data frame containing the information of the variables.

- horizons:

  A vector of time horizons which we want the absolute risk predictions
  to be evaluated at.

- event:

  The value for event of interest which we want the absolute risk
  predictions to be evaluated for. This must be one of the values in the
  status variable of the data.

- rhs:

  A right hand sided formula indicating the variables to be used in
  estimating the inverse probability of censoring weighting (IPCW)
  model. Default is `~1`.

- method:

  Resampling method to be used for hyper-parameter tuning. Values can
  be: `'cv'` for cross validation, `'repcv'` for repeated cross
  validation, `'lgocv'` for monte-carlo cross validation, `'loocv'` for
  leave one out cross validation and `'boot'` for bootstrap. Default is
  `'cv'`.

- k:

  Number of folds. Only applicable for `method='cv'` and
  `method='repcv'`. Default is 10.

- times:

  Repeat number of the resampling process. Only applicable for
  `method='repcv'`, `method='lgocv'` and `method='boot'`. Default is 25.

- p:

  The fraction of data to be used as the training set during resampling.
  Only applicable for `method='lgocv'`. Default is 0.7.

- strat.var:

  A single character indicating name of the strata variable to be used
  to create resamples. If numerical, groups will be specified based on
  percentiles. Default is `NULL` which considers status variable as a
  factor and creates the resamples based on different levels of it.

- metrics:

  Evaluation metric (loss function) to be used. Values can be `'Brier'`
  for IPCW brier score, `'AUC'` for IPCW AUC or a vector of both.
  Default is `'Brier'`.

- final.metric:

  The evaluation metric to decide the best hyper-parameters set for the
  final fits on the whole data. When `NULL` which is the default value,
  it takes the value from `metrics`. If both `'Brier'` and `'AUC'` were
  specified in metrics and `final.metric` is `NULL`, `'Brier'` will be
  used.

- alpha.grid:

  A named list containing a sequence of alpha values to be evaluated for
  each cause-specific model. Names of the list must be the events and
  exactly the same as values in the status variable. Default is `NULL`
  which orders the function to set `seq(0,1,.5)` for all cause-specific
  models. See \`Details\` for more information.

- lambda.grid:

  A named list containing a sequence of lambda values to be evaluated
  for each cause-specific model. Names of the list must be the events
  and exactly the same as values in the status variable. Default is
  `NULL` which orders the function to calculate exclusive lambda
  sequences for all causes. See \`Details\` for more information.

- nlambdas.list:

  A names list of single integers indicating the length of lambda
  sequences which are calculated automatically by the function for each
  cause. Only applicable when `lambda.grid=NULL`. Default is NULL which
  sets all lengths to 5. See \`Details\` for more information.

- grow.by:

  Difference between the values in the growing sequence of lambda values
  to find the maximum value that makes the null model. Only applicable
  when `lambda.grid=NULL`. Default is 0.01. See \`Details\` for more
  information.

- standardize:

  Logical indicating whether the variables must be standardized or not
  during model fitting procedures. Default is `TRUE`.

- keep:

  A character vector of the names of variables that should not be shrunk
  in all model fitting procedures. Default is `NULL`.

- preProc.fun:

  A function that accepts a data and returns a modified version of it
  that has gone through the user's desired pre-processing steps. All
  modifications from this function will be done during the resampling
  procedures to avoid data leakage. It will modify all training and test
  set(s) during the validation unless other argument `preProc.fun.test`
  is specified by user and then it only affects the training set(s).
  Default is `function(x) x`. Also see the description of
  `preProc.fun.test` argument.

- preProc.fun.test:

  A function the exact same characteristics and description as
  `preProc.fun` argument. If user specifies a separate function for
  `preProc.fun.test`, it will only affect test set(s) during validation
  while the function from `preProc.fun` will affect the training set(s).
  Default is `NULL` which means function from `preProc.fun` will be used
  on both training and test set(s) during validation. Also see the
  description of `preProc.fun` argument.

- parallel:

  Logical indicating whether the tuning process should be performed in
  parallel or not. Default is `FALSE`.

- preProc.pkgs:

  A character vector containing the names of packages that was used in
  creating user's `preProc.fun` while using parallel computation. Only
  applicable if `parallel=T` and `preProc.fun` is a user specified
  function using functions from other packages. See 'Examples' for
  details.

- preProc.globals:

  A character vector containing names of objects included in
  `preProc.fun` to be considered as global objects while using parallel
  computation. The most frequent ones are the names of the user
  specified pre processing function or functions within this function.
  Only applicable if `parallel=T` and `preProc.fun` is a user specified
  function. See 'Examples' for details.

- core.nums:

  Number of CPU cores to be used for parallel computation. Only
  applicable if `parallel=T`. Default is `future::availableCores()/2`.

## Value

A list containing the detailed information of the hyper-parameter tuning
and the validation process, best combination of hyper-parameters and the
final fits based on the whole data using the best obtained
hyper-parameters. Use `$` to explore all the involved information.

## Details

`tune_penCSC` has the ability to automatically determine the candidate
sequences of alpha & lambda values. Setting any of `alpha.grid` &
`lambda.grid` to `NULL` will order the function to calculate them
automatically. The process of determining the lambda values
automatically is by:

1.  Starting from lambda=0, the algorithm fits LASSO models until
    finding a lambda value that creates a NULL model where all variables
    were shrunk to be exactly zero.

2.  The obtained lambda value will be used as the maximum value of a
    sequence starting from 0. The length of this sequence is controlled
    by values in `nlambdas.list`.

This will be done for each cause-specific model to create exclusive
sequences of lambdas for each of them.

## References

Friedman J, Hastie T, Tibshirani R (2010). "Regularization Paths for
Generalized Linear Models via Coordinate Descent." Journal of
Statistical Software, 33(1), 1-22.
[doi:10.18637/jss.v033.i01](https://doi.org/10.18637/jss.v033.i01) ,
<https://www.jstatsoft.org/v33/i01/>.

Saadati, M, Beyersmann, J, Kopp-Schneider, A, Benner, A. Prediction
accuracy and variable selection for penalized cause-specific hazards
models. Biometrical Journal. 2018; 60: 288– 306.
[doi:10.1002/bimj.201600242](https://doi.org/10.1002/bimj.201600242) .

Gerds TA, Kattan MW (2021). Medical Risk Prediction Models: With Ties to
Machine Learning (1st ed.). Chapman and Hall/CRC.
[doi:10.1201/9781138384484](https://doi.org/10.1201/9781138384484)

Pfeiffer, R. M., & Gail, M. M. (2017). Absolute risk: Methods and
applications in clinical management and public health.

Kuhn, M. (2008). Building Predictive Models in R Using the caret
Package. Journal of Statistical Software, 28(5), 1–26.
[doi:10.18637/jss.v028.i05](https://doi.org/10.18637/jss.v028.i05) .

Bengtsson H (2021). “A Unifying Framework for Parallel and Distributed
Processing in R using Futures.” The R Journal, 13(2), 208–227.
[doi:10.32614/RJ-2021-048](https://doi.org/10.32614/RJ-2021-048) .

Vaughan D, Dancho M (2022). furrr: Apply Mapping Functions in Parallel
using Futures. <https://github.com/DavisVaughan/furrr>,
<https://furrr.futureverse.org/>.

Therneau T (2022). A Package for Survival Analysis in R. R package
version 3.3-1, <https://CRAN.R-project.org/package=survival>.

Wickham H, Averick M, Bryan J, Chang W, McGowan L, François R, et al.
Welcome to the tidyverse. J Open Source Softw. 2019 Nov 21;4(43):1686.

Bache S, Wickham H (2022). magrittr: A Forward-Pipe Operator for R.
<https://magrittr.tidyverse.org>,
<https://github.com/tidyverse/magrittr>.

## Author

Shahin Roshani

## Examples

``` r
# \donttest{

library(survival)

library(riskRegression)

data(Melanoma)

vl <- list('1'=~age+sex+epicel+ici,

          '2'=c('age','ulcer','thick','invasion'))

al <- list('1'=0,'2'=c(.5,1))

#External standardization function with data frame as its input and output

library(recipes)
#> 
#> Attaching package: ‘recipes’
#> The following object is masked from ‘package:stringr’:
#> 
#>     fixed
#> The following object is masked from ‘package:stats’:
#> 
#>     step

std.fun <- function(data){

 cont_vars <- data %>% select(where(~is.numeric(.))) %>% names

 cont_vars <- cont_vars[-which(cont_vars %in% c('time','status'))]

 #External functions from recipes package are being used

 recipe(~.,data=data) %>%

   step_center(all_of(cont_vars)) %>%

   step_scale(all_of(cont_vars)) %>%

   prep(training=data) %>% juice

}

set.seed(233)

test <- tune_penCSC(time='time',status='status',vars.list=vl,data=Melanoma,horizons=1825,

                   event=1,method='cv',k=5,metrics='AUC',alpha.grid=al,standardize=FALSE,

                   preProc.fun=std.fun,parallel=TRUE,preProc.pkgs='recipes')
#> 
#> Process was done in 18.12464 secs.

test
#> $`1825`
#> $`Event: 1`
#> 6 x 1 sparse Matrix of class "dgCMatrix"
#>                        1
#> age            0.3748226
#> sexMale        0.7967235
#> epicelpresent -1.1097122
#> ici1           1.7729566
#> ici2           1.8709727
#> ici3           2.5489947
#> 
#> $`Event: 2`
#> 5 x 1 sparse Matrix of class "dgCMatrix"
#>                         1
#> age             0.3122109
#> ulcerpresent    .        
#> thick           .        
#> invasionlevel.1 .        
#> invasionlevel.2 .        
#> 
#> 

# }
```
