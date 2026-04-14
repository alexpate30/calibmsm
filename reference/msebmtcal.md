# European Group for Blood and Marrow Transplantation data in `msdata` format.

The [`ebmt4`](https://rdrr.io/pkg/mstate/man/EBMT-data.html) data
converted into `msdata` format (see
[`msprep`](https://rdrr.io/pkg/mstate/man/msprep.html)), using the
processes implemented in the `mstate` package. Code for the derivation
of this dataset is provided in the source code for the package.

## Usage

``` r
msebmtcal
```

## Format

### 'msebmtcal'

A data frame in `msdata` format (see
[`msprep`](https://rdrr.io/pkg/mstate/man/msprep.html)) with 15,512 rows
and 8 columns:

- id:

  Patient indentifier

- from:

  transition from state

- to:

  transition to state

- trans:

  transition number

- Tstart:

  time entered state 'from'

- Tstop:

  time leaving state 'from'

- time:

  time in state 'from'

- status:

  event indicator, 1:transitioned to state 'to'

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
