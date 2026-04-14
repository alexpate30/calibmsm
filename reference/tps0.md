# Predicted transition probabilities out of transplant state made at time s = 0

Data frame containing the predicted transition probabilities out of
state j = 1 made at time s = 0. The predicted transition probabilities
were estimated by fitting a multistate model to the
[`ebmt4`](https://rdrr.io/pkg/mstate/man/EBMT-data.html) data using a
leave-one-out approach. Code for deriving this dataset is provided in
the source code for `calibmsm`. Code for the derivation of this dataset
is provided in the source code for the package.

## Usage

``` r
tps0
```

## Format

### 'tps0'

A data frame with 13,674 (CHANGE) rows and 14 columns:

- id:

  Patient indentifier

- pstate1, pstate2, pstate3, pstate4, pstate5, pstate6:

  Predicted transition probabilities of transitions into states 1 to 6

- se1, se2, se3, se4, se5, se6:

  Standard error of the predicted transition probabilities of
  transitions into states 1 to 6

- j:

  State from which the predicted transition probabilities are estimated
  from

## Source

This dataset was derived from data made available within the `mstate`
package, see [`ebmt4`](https://rdrr.io/pkg/mstate/man/EBMT-data.html).
The data was originally provided by the European Group for Blood and
Marrow Transplantation (https://www.ebmt.org/). We reiterate the source
statement given by the developers of `mstate`: "We acknowledge the
European Society for Blood and Marrow Transplantation (EBMT) for making
available these data. Disclaimer: these data were simplified for the
purpose of illustration of the analysis of competing risks and
multi-state models and do not reflect any real life situation. No
clinical conclusions should be drawn from these data."

## References

EBMT (2023). “Data from the European Society for Blood and Marrow
Transplantation.” URL
https://search.r-project.org/CRAN/refmans/mstate/html/EBMT-data.html.

de Wreede LC, Fiocco M, Putter H (2011). “mstate: An R Package for the
Analysis of Competing Risks and Multi-State Models.” *Journal of
Statistical Software*, 38(7).
