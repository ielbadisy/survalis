# Plot an interpretation

Draws the result of
[`interpret()`](https://elbadisyimad.com/survalis/reference/interpret.md)
with the plot of its method. Plots have no title unless `title` is
given, and return a ggplot to which layers can be added.

## Usage

``` r
# S3 method for class 'survalis_interpret'
plot(x, ...)
```

## Arguments

- x:

  A `survalis_interpret` object.

- ...:

  Passed to the method's plot. Useful ones:

  - `title` (all methods; `subtitle` and `caption` too for `"roc"` and
    `"dca"`).

  - `time.labels` (`"pdp"`, `"ice"`, `"ale"`): labels of the time panels
    or legend. `NULL` gives "t = 365", a string is the unit (`"days"`
    gives "t = 365 days"), and a function maps the times to labels, for
    example `function(t) paste(t / 365, "years")`.

  - `colors`: `"pdp"` / `"ice"` the ICE and PDP colours (default grey
    and blue); `"roc"` the curve colour; `"dca"` the model, treat-all
    and treat-none colours. Colours are set through the theme palette
    (ggplot2 \>= 4.0), so a `scale_colour_manual()` added afterwards
    also works, without a message; a complete theme added afterwards
    (e.g. `theme_bw()`) resets them, so prefer `colors` then.

  - `auc` (`"roc"`): print the AUC in the panel (default `TRUE`).

## Value

A ggplot.

## Examples

``` r
m <- fit(Surv(time, status) ~ age + karno, veteran, "coxph")
plot(interpret(m, "pdp", feature = "karno", times = c(90, 180)),
     time.labels = "days")

plot(interpret(m, "dca", eval_time = 180), colors = c("black", "grey50", "grey80"))
```
