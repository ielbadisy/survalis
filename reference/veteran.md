# Veteran's Administration Lung Cancer Trial Data

This is the `veteran` dataset originally from the survival package,
containing data from a randomized trial of lung cancer treatments.

## Usage

``` r
veteran
```

## Format

A data frame with 137 observations and 8 variables:

- trt:

  Treatment: 1=standard, 2=test

- celltype:

  Cell type: squamous, smallcell, adeno, large

- time:

  Survival time in days

- status:

  Censoring status: 1=dead, 0=alive

- karno:

  Karnofsky performance score (higher = better)

- diagtime:

  Months from diagnosis to randomization

- age:

  Age in years

- prior:

  Prior therapy: 0=no, 10=yes

## Source

survival package, originally from Kalbfleisch and Prentice (1980) *The
Statistical Analysis of Failure Time Data*.

## Examples

``` r
head(veteran)
#>   trt celltype time status karno diagtime age prior
#> 1   1 squamous   72      1    60        7  69     0
#> 2   1 squamous  411      1    70        5  64    10
#> 3   1 squamous  228      1    60        3  38     0
#> 4   1 squamous  126      1    60        9  63    10
#> 5   1 squamous  118      1    70       11  65    10
#> 6   1 squamous   10      1    20        5  49     0
summary(veteran$time)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>     1.0    25.0    80.0   121.6   144.0   999.0 
```
