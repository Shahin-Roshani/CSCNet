# predict.penCSC

Flexible prediction method for the objects of class \`penCSC\` including
the absolute risk prediction.

## Usage

``` r
# S3 method for class 'penCSC'
predict(object, newX, event = NULL, time, type = "lp", reference = "zero", ...)
```

## Arguments

- object:

  An object of class \`penCSC\`.

- newX:

  A data frame containing the information of variables related to new
  records. Information of variables not included in the model creation
  will be ignored.

- event:

  A vector of event codes which we want predictions for. This must be
  the same as values in the status variable of the data that was used to
  create the models. If `NULL`, absolute risk will be calculated for all
  involved events. Default is `NULL` which returns values for all
  involved causes.

- time:

  A vector of time horizons which we want absolute risk predictions at.
  Only applicable when `type='absRisk'`.

- type:

  Type of the predictions. Valid values are: `'lp'` or `'link'` for
  linear predictors, `'risk'` or `'response'` for `exp(lp)` and finally
  `'absRisk'` for semi-parametric estimates of absolute risk.

- reference:

  Reference for centering predictions. Valid values are `'zero'` and
  `'sample'`. Default is `'zero'`. For more information on referencing
  see details in `?predict.coxph`.

- ...:

  Additional arguments. Not used by `predict.penCSC`.

## Value

A tibble containing the predictions based on the input arguments.

## References

Pfeiffer, R. M., & Gail, M. M. (2017). Absolute risk: Methods and
applications in clinical management and public health.

Aalen, O.O. (1978) Nonparametric Inference for a Family of Counting
Processes. The Annals of Statistics, 6, 701-726.
[doi:10.1214/aos/1176344247](https://doi.org/10.1214/aos/1176344247) .

Wickham H, Averick M, Bryan J, Chang W, McGowan L, François R, et al.
Welcome to the tidyverse. J Open Source Softw. 2019 Nov 21;4(43):1686.

Bache S, Wickham H (2022). magrittr: A Forward-Pipe Operator for R.
<https://magrittr.tidyverse.org>,
<https://github.com/tidyverse/magrittr>.

Friedman J, Hastie T, Tibshirani R (2010). "Regularization Paths for
Generalized Linear Models via Coordinate Descent." Journal of
Statistical Software, 33(1), 1-22.
[doi:10.18637/jss.v033.i01](https://doi.org/10.18637/jss.v033.i01) ,
<https://www.jstatsoft.org/v33/i01/>.

## Author

Shahin Roshani

## Examples

``` r
library(riskRegression)

data(Melanoma)

vl <- list('1'=c('age','sex','ulcer','thick'),

          '2'=~age+sex+epicel+thick+ici)

al <- list('1'=0,'2'=.5)

ll <- list('1'=.01,'2'=.04)

penfit <- penCSC(time='time',status='status',vars.list=vl,

                data=Melanoma,alpha.list=al,lambda.list=ll)

predict(penfit,Melanoma[1:5,],type='lp')
#> # A tibble: 10 × 3
#>       id event prediction
#>    <int> <chr>      <dbl>
#>  1     1 1           3.17
#>  2     2 1           1.16
#>  3     3 1           1.05
#>  4     4 1           1.15
#>  5     5 1           3.46
#>  6     1 2           2.66
#>  7     2 2           1.96
#>  8     3 2           1.43
#>  9     4 2           2.48
#> 10     5 2           1.82

predict(penfit,Melanoma[1:5,],type='response')
#> # A tibble: 10 × 3
#>       id event prediction
#>    <int> <chr>      <dbl>
#>  1     1 1          23.9 
#>  2     2 1           3.18
#>  3     3 1           2.87
#>  4     4 1           3.17
#>  5     5 1          32.0 
#>  6     1 2          14.2 
#>  7     2 2           7.08
#>  8     3 2           4.19
#>  9     4 2          12.0 
#> 10     5 2           6.15

predict(penfit,Melanoma[1:5,],type='absRisk',event=1:2,time=1825*(1:2))
#> # A tibble: 20 × 4
#>       id event horizon absoluteRisk
#>    <int> <int>   <dbl>        <dbl>
#>  1     1     1    1825       0.599 
#>  2     2     1    1825       0.125 
#>  3     3     1    1825       0.115 
#>  4     4     1    1825       0.121 
#>  5     5     1    1825       0.720 
#>  6     1     2    1825       0.0725
#>  7     2     2    1825       0.0451
#>  8     3     2    1825       0.0271
#>  9     4     2    1825       0.0749
#> 10     5     2    1825       0.0302
#> 11     1     1    3650       0.763 
#> 12     2     1    3650       0.207 
#> 13     3     1    3650       0.193 
#> 14     4     1    3650       0.200 
#> 15     5     1    3650       0.868 
#> 16     1     2    3650       0.114 
#> 17     2     2    3650       0.140 
#> 18     3     2    3650       0.0879
#> 19     4     2    3650       0.222 
#> 20     5     2    3650       0.0423
```
