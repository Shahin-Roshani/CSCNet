# CSCNet 0.1.1

In `predict.penCSC`:

* Corrected the issue in setting factor levels to character variables in `newdata` argument.

In `tune_penCSC`:

* Added `strat.var` argument to specify the variable to be used as the grouping variable when creating the resamples. With this, the issue that numerical status variable was not considered as a factor to create groups was also fixed.

* Added `preProc.fun.test` argument to specify a separate pre-processing function to affect test set(s) during validation. This allows the user to make a distinction between pre-processing steps of training and test set(s) if needed.

Also, `as.character.POSIXt()` was replaced with `format()` in different parts following the latest change in the behavior of the former function.
