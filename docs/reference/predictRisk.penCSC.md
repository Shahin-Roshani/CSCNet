# predictRisk.penCSC

predictRisk method for absolute risk prediction. This is mainly for
compatibility of 'CSCNet' with functions of 'riskRegression' package.

## Usage

``` r
# S3 method for class 'penCSC'
predictRisk(object, newdata, times, cause, ...)
```

## Arguments

- object:

  An object of class 'penCSC'.

- newdata:

  A data frame containing the variable information of new records.

- times:

  A vector of time horizons which we want the absolute risk predictions
  at.

- cause:

  A single value indicating the event of interest which we want the
  absolute risk predictions for. This value should be one of the values
  in the status variable of the data.

- ...:

  Additional arguments. Not used by `predictRisk.penCSC`.

## Value

A matrix with columns of absolute risk predictions of individuals for
each requested time horizon.

## References

Wickham H, Averick M, Bryan J, Chang W, McGowan L, François R, et al.
Welcome to the tidyverse. J Open Source Softw. 2019 Nov 21;4(43):1686.

Bache S, Wickham H (2022). magrittr: A Forward-Pipe Operator for R.
<https://magrittr.tidyverse.org>,
<https://github.com/tidyverse/magrittr>.

## See also

<https://www.rdocumentation.org/packages/riskRegression/versions/1.3.7/topics/predictRisk>

Details in: <https://rdrr.io/cran/riskRegression/man/Score.html>

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

predictRisk(penfit,Melanoma[1:5,],times=1825*(1:2),cause=1)
#>           1825      3650
#> [1,] 0.5989514 0.7634726
#> [2,] 0.1245976 0.2072521
#> [3,] 0.1147573 0.1931740
#> [4,] 0.1213908 0.1996565
#> [5,] 0.7199651 0.8683818
```
