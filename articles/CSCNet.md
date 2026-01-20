# CSCNet vignette

CSCNet is package with flexible tools for fitting and evaluating
cause-specific cox models with elastic-net penalty. Each cause is
modeled in a separate penalized cox model (using elastic-net penalty)
with its exclusive $\alpha$ and $\lambda$ assuming other involved
competing causes as censored.

### Regularized cause-specific cox and absolute risk predictions

In this package we will use a simulated data using ‘riskRegression’
package (which will load up with ‘CSCNet’).

``` r
library(CSCNet)

library(riskRegression)

set.seed(123)

d <- sampleData(n = 500,outcome = 'competing.risks')

knitr::kable(head(d),digits=3)
```

|     X6 |     X7 |     X8 |     X9 |    X10 | X1  | X2  | X3  | X4  | X5  | eventtime1 | eventtime2 | censtime |   time | event |
|-------:|-------:|-------:|-------:|-------:|:----|:----|:----|:----|:----|-----------:|-----------:|---------:|-------:|------:|
| 38.651 | 62.599 | -1.117 | -1.253 | -0.702 | 0   | 1   | 0   | 1   | 0   |     17.595 |     10.196 |   11.164 | 10.196 |     2 |
| 75.335 | 59.612 | -0.892 | -0.111 |  0.882 | 1   | 0   | 0   | 1   | 1   |      1.652 |     10.020 |    4.877 |  1.652 |     1 |
| 70.317 | 64.125 |  0.872 | -1.413 | -0.133 | 0   | 0   | 0   | 0   | 0   |      2.006 |     11.171 |    9.456 |  2.006 |     1 |
| 55.388 | 65.603 |  1.869 | -1.983 | -1.121 | 0   | 0   | 0   | 0   | 1   |      2.082 |     15.894 |    3.527 |  2.082 |     1 |
| 59.704 | 59.389 | -0.124 |  0.784 |  0.461 | 0   | 0   | 0   | 1   | 1   |      8.961 |     10.033 |    2.477 |  2.477 |     0 |
| 67.326 | 69.727 |  0.107 |  0.901 |  1.524 | 0   | 0   | 0   | 0   | 0   |     23.113 |     13.137 |   17.575 | 13.137 |     2 |

``` r
table(d$event)

  0   1   2 
140 225 135 
```

There are 2 events in the data coded as 1 & 2. To introduce how setting
up variables and hyper-parameters works in CSCNet, we will fit the a
model with the following hyper-parameters:
$$\left( \alpha_{1},\alpha_{2},\lambda_{1},\lambda_{2} \right) = (0,0.5,0.01,0.02)$$
We set variables affecting the event: 1 as `X1, X3, X7, X9, X10` and
variables affecting event: 2 as `X1, X2, X6, X10`.

#### Fitting regularized cause-specific cox models

In CSCNet, setting variables and hyper-parameters are done through named
lists. Variables and hyper-parameters related to each involved cause are
stored in list positions with the name of that position being that
cause. Of course these names must be the same as values in the status
variable in the data.

``` r
vl <- list('1'=~X1+X3+X7+X9+X10,

           '2'=c('X1','X2','X6','X10'))

penfit <- penCSC(time = 'time',status = 'event',vars.list = vl,data = d,
                 
                 alpha.list = list('1'=0,'2'=.5),lambda.list = list('1'=.01,'2'=.02))

penfit
$`Event: 1`
5 x 1 sparse Matrix of class "dgCMatrix"
              1
X11  0.97724748
X31  0.17673877
X7  -0.04914618
X9  -0.52737860
X10 -0.09007295

$`Event: 2`
4 x 1 sparse Matrix of class "dgCMatrix"
              1
X11 -0.71082304
X21  0.28319963
X6   .         
X10  0.09318867
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
predict(penfit,d[1:3,],type='lp',event=1) %>% as.data.frame
  id event prediction
1  1     1  -2.352336
2  2     1  -1.973192
3  3     1  -2.394408
```

Or the risk values of the same individuals for all involved causes:

``` r
predict(penfit,d[1:3,],type='response') %>% as.data.frame
  id event prediction
1  1     1 0.09514665
2  2     1 0.13901236
3  3     1 0.09122671
4  1     2 1.24337231
5  2     2 0.53333328
6  3     2 0.98764831
```

Now let’s say we want estimates of absolute risks related to the event:
1 as our event of interest at median and 3rd quartile of the follow-up
times:

``` r
predict(penfit,d[1:3,],type='absRisk',event=1,time=summary(d$time)[c(3,5)]) %>% as.data.frame
  id event  horizon absoluteRisk
1  1     1 3.846404    0.4212294
2  2     1 3.846404    0.5625537
3  3     1 3.846404    0.4124798
4  1     1 6.243125    0.5451776
5  2     1 6.243125    0.7129566
6  3     1 6.243125    0.5423173
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

1.  The algorithm fits LASSO models until finding a $\lambda$ value that
    creates a NULL model where all variables were shrunk to be
    exactly 0. The path to reach maximum $\lambda$ is specified
    internally by `glmnet`.
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
Let’s say we want to tune our model for absolute risk prediction of
event: 1 at the median follow-up time based on time dependent (IPCW) AUC
as the loss function (evaluation metric) through a 5-fold cross
validation process:

**Note:** The following code chunk is not executed during vignette
building, as its output depends on random number generation that may
vary across platforms. Users are encouraged to run the code locally and
explore different random seeds.

``` r
#Function to standardize numerical predictors using functions from recipes package

library(recipes)

pp.fun <- function(data){

  recipe(time+event~.,data=data) %>% 
    
    step_center(all_numeric_predictors()) %>% 
    
    step_scale(all_numeric_predictors()) %>% 
    
    prep(training=data) %>% 
    
    bake(new_data=NULL)

}

set.seed(123)

tri.l <- caret::createFolds(as.factor(d$event),k=5,list=T,returnTrain=T)

tune_obj <- tune_penCSC(time = 'time',status = 'event',vars.list = vl,data = d,
                        
                        horizons = median(d$time),event = 1,tri.list = tri.l,metrics = 'AUC',
                        
                        alpha.grid = list('1'=0,'2'=c(.5,1)),preProc.fun = pp.fun,
                        
                        standardize = F,parallel = T,preProc.pkgs = 'recipes')

tune_obj$validation_result %>% arrange(desc(mean.AUC)) %>% head

tune_obj$final_params

tune_obj$final_fits
```
