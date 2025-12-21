# penCSC

Function to fit penalized cause-specific-cox with elastic-net penalty.

## Usage

``` r
penCSC(
  time,
  status,
  vars.list,
  data,
  alpha.list,
  lambda.list,
  standardize = TRUE,
  keep = NULL
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

- alpha.list:

  A named list containing the single alpha values of each cause-specific
  model. Names of the list must be the events and exactly the same as
  values in the status variable. See \`Examples\` for details.

- lambda.list:

  A named list containing the single lambda values of each
  cause-specific model. Names of the list must be the events and exactly
  the same as values in the status variable. See \`Examples\` for
  details.

- standardize:

  Logical indicating whether the variables must be standardized or not.
  Default is `TRUE`.

- keep:

  A character vector of the names of variables that should not be
  shrunk. Default is `NULL`.

## Value

A named list containing all the information related to the used data and
the fitted models for all causes. Use `$` to explore all the involved
information.

## References

Friedman J, Hastie T, Tibshirani R (2010). "Regularization Paths for
Generalized Linear Models via Coordinate Descent." Journal of
Statistical Software, 33(1), 1-22.
[doi:10.18637/jss.v033.i01](https://doi.org/10.18637/jss.v033.i01) ,
<https://www.jstatsoft.org/v33/i01/>.

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
library(riskRegression)
#> riskRegression version 2025.09.17

data(Melanoma)

vl <- list('1'=c('age','sex','ulcer','thick'),

          '2'=~age+sex+epicel+thick+ici)

al <- list('1'=0,'2'=.5)

ll <- list('1'=.01,'2'=.04)

penCSC(time='time',status='status',vars.list=vl,

      data=Melanoma,alpha.list=al,lambda.list=ll)
#> $`Event: 1`
#> 4 x 1 sparse Matrix of class "dgCMatrix"
#>                       1
#> age          0.01184358
#> sexMale      0.42366136
#> ulcerpresent 1.11871164
#> thick        0.10816182
#> 
#> $`Event: 2`
#> 7 x 1 sparse Matrix of class "dgCMatrix"
#>                        1
#> age           0.03494673
#> sexMale       .         
#> epicelpresent .         
#> thick         .         
#> ici1          .         
#> ici2          .         
#> ici3          .         
#> 
```
