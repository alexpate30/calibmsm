# European Group for Blood and Marrow Transplantation data (one row per individual)

Used in vignette/article:
Comparison-with-graphical-calibration-curves-in-competing-risks-setting.

## Usage

``` r
ebmtcal_cmprsk
```

## Format

### 'ebmtcal_cmprsk'

A data frame with 2,279 rows and 17 columns:

- id:

  Patient indentifier

- rec, rec.s:

  Time until and event indicator for recovery variable

- ae, ae.s:

  Time until and event indicator for adverse event variable

- recae, recae.s:

  Time until and event indicator for recovery + adverse event variable

- rel, rel.s:

  Time until and event indicator for relapse variable

- srv, srv.s:

  Time until and event indicator for death variable

- year:

  Year of transplant

- agecl:

  Age at transplant

- proph:

  Prophylaxis

- match:

  Donor-recipient match

- dtcens:

  Time of censoring

- dtcens_s:

  Event indicator, 1:censoring occured, 0: absorbing state entered
  before censoring occured

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

## Details

A data frame of 2,279 individuals with blood cancer who have undergone a
transplant. This data is identical to the
[`ebmt4`](https://rdrr.io/pkg/mstate/man/EBMT-data.html) data, except
two extra variables have been derived, time until censoring and a
censoring indicator, which are required to assess calibration using some
of the methods in `calibmsm`. Specifically, the time until censoring ar
calculated in the setting of a competing risks model out of the first
state, where no further transitions can be made. This means entry into
any state (as they are all absorbing states) will have the effect of
preventing censoring from occurring, and `dtcens` and `dtcens_s` will be
different than the values found in
[`ebmtcal`](https://alexpate30.github.io/calibmsm/reference/ebmtcal.md).
This dataset has been designed to be used alongside dataset
[`msebmtcal_cmprsk`](https://alexpate30.github.io/calibmsm/reference/msebmtcal_cmprsk.md),
when assessing the calibration of a competing risks model. Code for the
derivation of this dataset is provided in the source code for the
package.

## References

EBMT (2023). “Data from the European Society for Blood and Marrow
Transplantation.” URL
https://search.r-project.org/CRAN/refmans/mstate/html/EBMT-data.html.

de Wreede LC, Fiocco M, Putter H (2011). “mstate: An R Package for the
Analysis of Competing Risks and Multi-State Models.” *Journal of
Statistical Software*, 38(7).
