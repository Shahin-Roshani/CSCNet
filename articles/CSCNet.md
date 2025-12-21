# CSCNet vignette

CSCNet is package with flexible tools for fitting and evaluating
cause-specific cox models with elastic-net penalty. Each cause is
modeled in a separate penalized cox model (using elastic-net penalty)
with its exclusive $\alpha$ and $\lambda$ assuming other involved
competing causes as censored.

### Regularized cause-specific cox and absolute risk predictions

In this package we will use `Melanoma` data from ‘riskRegression’
package (which will load up with ‘CSCNet’) so we start by loading the
package and the `Melanoma` data.

``` r
library(CSCNet)
library(riskRegression)
data(Melanoma)
as_tibble(Melanoma)
 [38;5;246m# A tibble: 205 × 11 [39m
    time status event     invasion ici   epicel ulcer thick sex     age logthick
    [3m [38;5;246m<int> [39m [23m   [3m [38;5;246m<dbl> [39m [23m  [3m [38;5;246m<fct> [39m [23m      [3m [38;5;246m<fct> [39m [23m     [3m [38;5;246m<fct> [39m [23m  [3m [38;5;246m<fct> [39m [23m   [3m [38;5;246m<fct> [39m [23m  [3m [38;5;246m<dbl> [39m [23m  [3m [38;5;246m<fct> [39m [23m  [3m [38;5;246m<int> [39m [23m     [3m [38;5;246m<dbl> [39m [23m
 [38;5;250m 1 [39m    10      2 death.ot… level.1  2     prese… pres…  6.76 Male     76    1.91 
 [38;5;250m 2 [39m    30      2 death.ot… level.0  0     not p… not …  0.65 Male     56   - [31m0 [39m [31m. [39m [31m431 [39m
 [38;5;250m 3 [39m    35      0 censored  level.1  2     not p… not …  1.34 Male     41    0.293
 [38;5;250m 4 [39m    99      2 death.ot… level.0  2     not p… not …  2.9  Fema…    71    1.06 
 [38;5;250m 5 [39m   185      1 death.ma… level.2  2     prese… pres… 12.1  Male     52    2.49 
 [38;5;250m 6 [39m   204      1 death.ma… level.2  2     not p… pres…  4.84 Male     28    1.58 
 [38;5;250m 7 [39m   210      1 death.ma… level.2  2     prese… pres…  5.16 Male     77    1.64 
 [38;5;250m 8 [39m   232      1 death.ma… level.2  2     not p… pres… 12.9  Male     49    2.56 
 [38;5;250m 9 [39m   232      2 death.ot… level.1  3     not p… pres…  3.22 Fema…    60    1.17 
 [38;5;250m10 [39m   279      1 death.ma… level.0  2     not p… pres…  7.41 Fema…    68    2.00 
 [38;5;246m# ℹ 195 more rows [39m
table(Melanoma$status)

  0   1   2 
134  57  14 
```

There are 2 events in the Melanoma data coded as 1 & 2. To introduce how
setting up variables and hyper-parameters works in CSCNet, we will fit
the a model with the following hyper-parameters to the `Melanoma` data:
$$\left( \alpha_{1},\alpha_{2},\lambda_{1},\lambda_{2} \right) = (0,0.5,0.01,0.02)$$
We set variables affecting the event: 1 as `age,sex,invasion,thick` and
variables affecting event: 2 as `age,sex,epicel,ici,thick`.

#### Fitting regularized cause-specific cox models

In CSCNet, setting variables and hyper-parameters are done through named
lists. Variables and hyper-parameters related to each involved cause are
stored in list positions with the name of that position being that
cause. Of course these names must be the same as values in the status
variable in the data.

``` r
vl <- list('1'=c('age','sex','invasion','thick'),
     
     '2'=~age+sex+epicel+ici+thick)

penfit <- penCSC(time = 'time',
                 
                 status = 'status',
                 
                 vars.list = vl,
                 
                 data = Melanoma,
                 
                 alpha.list = list('1'=0,'2'=.5),
                 
                 lambda.list = list('1'=.01,'2'=.02))

penfit
$`Event: 1`
5 x 1 sparse Matrix of class "dgCMatrix"
                          1
age             0.008018578
sexMale         0.547580959
invasionlevel.1 0.756922406
invasionlevel.2 0.591044240
thick           0.118568171

$`Event: 2`
7 x 1 sparse Matrix of class "dgCMatrix"
                        1
age            0.04839997
sexMale        0.11419057
epicelpresent  0.16891622
ici1          -0.13501846
ici2           .         
ici3           .         
thick          0.03242932
```

`penfit` is a comprehensive list with all information related to the
data and fitted models in detail that user can access.

**Note:** As we saw, variable specification in `vars.list` is possible
in 2 ways which are introducing a vector of variable names or a one hand
sided formula for different causes.

#### Predictions and semi-parametric estimates of absolute risk

Now to obtain predictions, specially estimates of the absolute risks,
`predict.penCSC` method was developed so user can obtain different forms
of values in the easiest way possible. By this method on objects of
class `penCSCS` and for different involved causes, user can obtain
values for linear predictors (`type='lp'` or `type='link'`), exponential
of linear predictors (`type='risk'` or `type='response'`) and finally
semi-parametric estimates of absolute risks (`type='absRisk'`) at
desired time horizons.

**Note:** Default value for `event` argument in `predict.penCSC` is
`NULL`. If user leaves it as that, values for all involved causes will
be returned.

Values of linear predictors for event: 1 related to 1st five individuals
of the data:

``` r
predict(penfit,Melanoma[1:5,],type='lp',event=1)
 [38;5;246m# A tibble: 5 × 3 [39m
     id event prediction
   [3m [38;5;246m<int> [39m [23m  [3m [38;5;246m<chr> [39m [23m       [3m [38;5;246m<dbl> [39m [23m
 [38;5;250m1 [39m     1 1          2.72 
 [38;5;250m2 [39m     2 1          1.07 
 [38;5;250m3 [39m     3 1          1.79 
 [38;5;250m4 [39m     4 1          0.913
 [38;5;250m5 [39m     5 1          2.99 
```

Or the risk values of the same individuals for all involved causes:

``` r
predict(penfit,Melanoma[1:5,],type='response')
 [38;5;246m# A tibble: 10 × 3 [39m
      id event prediction
    [3m [38;5;246m<int> [39m [23m  [3m [38;5;246m<chr> [39m [23m       [3m [38;5;246m<dbl> [39m [23m
 [38;5;250m 1 [39m     1 1          15.1 
 [38;5;250m 2 [39m     2 1           2.93
 [38;5;250m 3 [39m     3 1           6.00
 [38;5;250m 4 [39m     4 1           2.49
 [38;5;250m 5 [39m     5 1          19.8 
 [38;5;250m 6 [39m     1 2          65.4 
 [38;5;250m 7 [39m     2 2          17.2 
 [38;5;250m 8 [39m     3 2           8.52
 [38;5;250m 9 [39m     4 2          34.1 
 [38;5;250m10 [39m     5 2          24.3 
```

Now let’s say we want estimates of absolute risks related to the event:
1 as our event of interest at 3 and 5 year time horizons:

``` r
predict(penfit,Melanoma[1:5,],type='absRisk',event=1,time=365*c(3,5))
 [38;5;246m# A tibble: 10 × 4 [39m
      id event horizon absoluteRisk
    [3m [38;5;246m<int> [39m [23m  [3m [38;5;246m<dbl> [39m [23m    [3m [38;5;246m<dbl> [39m [23m         [3m [38;5;246m<dbl> [39m [23m
 [38;5;250m 1 [39m     1     1     [4m1 [24m095       0.374 
 [38;5;250m 2 [39m     2     1     [4m1 [24m095       0.095 [4m2 [24m
 [38;5;250m 3 [39m     3     1     [4m1 [24m095       0.187 
 [38;5;250m 4 [39m     4     1     [4m1 [24m095       0.079 [4m8 [24m
 [38;5;250m 5 [39m     5     1     [4m1 [24m095       0.480 
 [38;5;250m 6 [39m     1     1     [4m1 [24m825       0.525 
 [38;5;250m 7 [39m     2     1     [4m1 [24m825       0.153 
 [38;5;250m 8 [39m     3     1     [4m1 [24m825       0.292 
 [38;5;250m 9 [39m     4     1     [4m1 [24m825       0.128 
 [38;5;250m10 [39m     5     1     [4m1 [24m825       0.654 
```

**Note:** There’s also `predictRisk.penCSC` to obtain absolute risk
predictions. This method was developed for compatibility with tools from
‘riskRegression’ package.

### Tuning the hyper-parameters

The above example was for illustration purposes. In real world analysis,
one must tune the hyper-parameters with respect to a proper loss
function through resampling procedures. `tune_penCSC` is a comprehensive
function that was built for this purpose on regularized cause-specific
cox models.

Like before, specification of variables and hyper-parameters are done
through named lists and sequences of candidate hyper-parameters related
to each involved cause are stored in list positions with the name of
that position being that cause. After that, `tune_penCSC` will create
all possible combinations from user’s specified sequences and evaluates
them using either IPCW brier score or IPCW AUC (as loss functions) based
on absolute risk predictions of the event of interest (linking) through
a chosen resampling process. Supported resampling procedures are: cross
validation (`method='cv'`), repeated cross validation
(`method='repcv'`), bootstrap (`method='boot'`), Monte-Carlo or leave
group out cross validation (`method='lgocv'`) and leave one out cross
validation (`method='loocv'`).

#### Automatic specification of hyper-parameters sequences

`tune_penCSC` has the ability to automatically determine the candidate
sequences of $\alpha$ & $\lambda$ values. Setting any of `alpha.grid` &
`lambda.grid` to `NULL` will order the function to calculate them
automatically.

While the automatic sequence of $\alpha$ values for all causes is
`seq(0,1,.5)`, the process of determining the $\lambda$ values
automatically is by:

1.  Starting from $\lambda = 0$, the algorithm fits LASSO models until
    finding a $\lambda$ value that creates a NULL model where all
    variables were shrunk to be exactly 0.
2.  The obtained $\lambda$ value will be used as the maximum value of a
    sequence starting from 0. The length of this sequence is controlled
    by values in `nlambdas.list`.

This will be done for each cause-specific model to create exclusive
sequences of $\lambda$s for each of them.

#### Pre-processing within resampling

If the data requires pre-processing steps, it must be done within the
resampling process to avoid data leakage. This can be achieved by using
`preProc.fun` argument of `tune_penCSC` function. This arguments accepts
a function that has a data as its only input and returns a modified
version of that data. Any pre-processing steps can be specified within
this function.

**Note:** `tune_penCSC` has the parallel processing option. If a user
has specified a function for pre-processing steps with global objects or
calls from other packages and wants to run the code in parallel, the
names of those extra packages and global objects must be given through
`preProc.pkgs` and `preProc.globals`.

Now let’s see all that was mentioned in this section in an example.
Let’s say we want to tune our model for 5 year absolute risk prediction
of event: 1 based on time dependent (IPCW) AUC as the loss function
(evaluation metric) through a 5-fold cross validation process:

``` r
#Writing a hypothetical pre-processing function

library(recipes)

Attaching package: 'recipes'
The following object is masked from 'package:stringr':

    fixed
The following object is masked from 'package:stats':

    step

std.fun <- function(data){

  cont_vars <- data %>% select(where(~is.numeric(.))) %>% names

  cont_vars <- cont_vars[-which(cont_vars %in% c('time','status'))]

  #External functions from recipes package are being used

  recipe(~.,data=data) %>%

    step_center(all_of(cont_vars)) %>%

    step_scale(all_of(cont_vars)) %>%

    prep(training=data) %>% juice

}

#Tuning a regularized cause-specific cox 

set.seed(455) #for reproducibility

tune_melanoma <- tune_penCSC(time = 'time',
                             
                             status = 'status',
                             
                             vars.list = vl,
                             
                             data = Melanoma,
                             
                             horizons = 365*5,
                             
                             event = 1,
                             
                             method = 'cv',
                             
                             k = 5,
                             
                             standardize = FALSE,
                             
                             metrics = 'AUC',
                             
                             alpha.grid = list('1'=0,'2'=c(.5,1)),
                             
                             preProc.fun = std.fun,
                             
                             parallel = TRUE,
                             
                             preProc.pkgs = 'recipes')

Process was done in 41.01057 secs.

tune_melanoma$validation_result %>% arrange(desc(mean.AUC)) %>% head
  alpha_1 alpha_2 lambda_1 lambda_2 horizon  mean.AUC
1       0     0.5   0.0425   0.0350    1825 0.7613930
2       0     0.5   0.1275   0.0525    1825 0.7336138
3       0     1.0   0.1275   0.0350    1825 0.7335694
4       0     1.0   0.1700   0.0700    1825 0.7324348
5       0     0.5   0.1700   0.0700    1825 0.7304677
6       0     1.0   0.1275   0.0700    1825 0.7280321

tune_melanoma$final_params
$`1825`
  alpha_1 alpha_2 lambda_1 lambda_2 horizon mean.AUC
1       0     0.5   0.0425    0.035    1825 0.761393

tune_melanoma$final_fits
$`1825`
$`Event: 1`
5 x 1 sparse Matrix of class "dgCMatrix"
                        1
age             0.1495339
sexMale         0.3396344
invasionlevel.1 0.3768998
invasionlevel.2 0.1294018
thick           0.4044303

$`Event: 2`
7 x 1 sparse Matrix of class "dgCMatrix"
                       1
age           0.63657394
sexMale       .         
epicelpresent .         
ici1          .         
ici2          .         
ici3          .         
thick         0.02163226
```
