# Changelog

## CSCNet 0.1.3

In `tune_penCSC`:

- Fixed an issue where the function only worked when the survival
  outcome included censoring.

## CSCNet 0.1.2

CRAN release: 2022-11-08

In `penCSC`:

- Added call to the results list.

In `predict.penCSC`:

- Corrected the issue that the `newdata` had to include all variables
  from the `data` when setting matching factor levels to characters in
  `newdata`.

In `tune_penCSC`:

- Corrected the issue that `keep` & `standardize` arguments were not
  affecting the `$final_fits` in the result.
