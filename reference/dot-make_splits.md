# Expand a resampling specification into concrete row-index splits

Expand a resampling specification into concrete row-index splits

## Usage

``` r
.make_splits(resampling, data, status_col = NULL)
```

## Arguments

- resampling:

  A `survalis_resampling` object.

- data:

  The (already row-filtered) modelling data frame.

- status_col:

  Name of the event-indicator column, used when
  `params$strata == "status"`.

## Value

A list of `list(analysis = <int>, assessment = <int>, id = <chr>)`.
