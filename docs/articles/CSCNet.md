# CSCNet vignette

CSCNet is package with flexible tools for fitting and evaluating
cause-specific cox models with elastic-net penalty. Each cause is
modeled in a separate penalized cox model (using elastic-net penalty)
with its exclusive $`\alpha`$ and $`\lambda`$ assuming other involved
competing causes as censored.

### Regularized cause-specific cox and absolute risk predictions

In this package we will use `Melanoma` data from ‘riskRegression’
package (which will load up with ‘CSCNet’) so we start by loading the
package and the `Melanoma` data.

``` r
library(CSCNet)
library(riskRegression)
data(Melanoma)
knitr::kable(head(Melanoma),digits=4)
```

| time | status | event | invasion | ici | epicel | ulcer | thick | sex | age | logthick |
|---:|---:|:---|:---|:---|:---|:---|---:|:---|---:|---:|
| 10 | 2 | death.other.causes | level.1 | 2 | present | present | 6.76 | Male | 76 | 1.9110 |
| 30 | 2 | death.other.causes | level.0 | 0 | not present | not present | 0.65 | Male | 56 | -0.4308 |
| 35 | 0 | censored | level.1 | 2 | not present | not present | 1.34 | Male | 41 | 0.2927 |
| 99 | 2 | death.other.causes | level.0 | 2 | not present | not present | 2.90 | Female | 71 | 1.0647 |
| 185 | 1 | death.malignant.melanoma | level.2 | 2 | present | present | 12.08 | Male | 52 | 2.4916 |
| 204 | 1 | death.malignant.melanoma | level.2 | 2 | not present | present | 4.84 | Male | 28 | 1.5769 |

``` r
table(Melanoma$status)

  0   1   2 
134  57  14 
```

There are 2 events in the Melanoma data coded as 1 & 2. To introduce how
setting up variables and hyper-parameters works in CSCNet, we will fit
the a model with the following hyper-parameters to the `Melanoma` data:
``` math
(\alpha_{1},\alpha_{2},\lambda_{1},\lambda_{2})=(0,0.5,0.01,0.02)
```
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

penfit <- penCSC(time = 'time',status = 'status',vars.list = vl,data = Melanoma,
                 
                 alpha.list = list('1'=0,'2'=.5),lambda.list = list('1'=.01,'2'=.02))

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

Values of linear predictors for event: 1 related to 1st three
individuals of the data:

``` r
predict(penfit,Melanoma[1:3,],type='lp',event=1) %>% as.data.frame
  id event prediction
1  1     1   2.715436
2  2     1   1.073691
3  3     1   1.792146
```

Or the risk values of the same individuals for all involved causes:

``` r
predict(penfit,Melanoma[1:3,],type='response') %>% as.data.frame
  id event prediction
1  1     1  15.111199
2  2     1   2.926159
3  3     1   6.002322
4  1     2  65.413371
5  2     2  17.213052
6  3     2   8.516833
```

Now let’s say we want estimates of absolute risks related to the event:
1 as our event of interest at 3 and 5 year time horizons:

``` r
predict(penfit,Melanoma[1:3,],type='absRisk',event=1,time=365*c(3,5)) %>% as.data.frame
  id event horizon absoluteRisk
1  1     1    1095   0.37363641
2  2     1    1095   0.09524797
3  3     1    1095   0.18730858
4  1     1    1825   0.52534831
5  2     1    1825   0.15302632
6  3     1    1825   0.29161813
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
sequences of $`\alpha`$ & $`\lambda`$ values. Setting any of
`alpha.grid` & `lambda.grid` to `NULL` will order the function to
calculate them automatically.

While the automatic sequence of $`\alpha`$ values for all causes is
`seq(0,1,.5)`, the process of determining the $`\lambda`$ values
automatically is by:

1.  The algorithm fits LASSO models until finding a $`\lambda`$ value
    that creates a NULL model where all variables were shrunk to be
    exactly 0. The path to reach maximum $`\lambda`$ is specified
    internally by `glmnet`.
2.  The obtained $`\lambda`$ value will be used as the maximum value of
    a sequence starting from 0. The length of this sequence is
    controlled by values in `nlambdas.list`.

This will be done for each cause-specific model to create exclusive
sequences of $`\lambda`$s for each of them.

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
#Function to standardize numerical predictors using functions from recipes package

library(recipes)

pp.fun <- function(data){

  recipe(time+status~.,data=data) %>% 
    
    step_center(all_numeric_predictors()) %>% 
    
    step_scale(all_numeric_predictors()) %>% 
    
    prep(training=data) %>% 
    
    bake(new_data=NULL)

}

set.seed(1331)

tri.l <- caret::createFolds(as.factor(Melanoma$status),k=3,list=T,returnTrain=T)

tune_melanoma <- tune_penCSC(time = 'time',status = 'status',vars.list = vl,data = Melanoma,
                             
                             horizons = 1095,event = 1,tri.list = tri.l,metrics = 'AUC',
                             
                             alpha.grid = list('1'=0,'2'=c(.5,1)),preProc.fun = pp.fun,
                             
                             standardize = F,parallel = T,preProc.pkgs = 'recipes')

tune_melanoma$validation_result %>% arrange(desc(mean.AUC)) %>% head
  alpha_1 alpha_2   lambda_1   lambda_2 horizon  mean.AUC
1       0     0.5 0.12668029 0.06461658    1095 0.7417364
2       0     1.0 0.12668029 0.04846244    1095 0.7411191
3       0     1.0 0.12668029 0.03230829    1095 0.7407155
4       0     0.5 0.12668029 0.03230829    1095 0.7400982
5       0     0.5 0.12668029 0.04846244    1095 0.7397054
6       0     1.0 0.08445352 0.03230829    1095 0.7396924

tune_melanoma$final_params
$`1095`
  alpha_1 alpha_2  lambda_1   lambda_2 horizon  mean.AUC
1       0     0.5 0.1266803 0.06461658    1095 0.7417364

tune_melanoma$final_fits
$`1095`
$`Event: 1`
5 x 1 sparse Matrix of class "dgCMatrix"
                        1
age             0.1417696
sexMale         0.1969270
invasionlevel.1 0.2063510
invasionlevel.2 0.0654537
thick           0.3556714

$`Event: 2`
7 x 1 sparse Matrix of class "dgCMatrix"
                      1
age           0.3523809
sexMale       .        
epicelpresent .        
ici1          .        
ici2          .        
ici3          .        
thick         .        
```
