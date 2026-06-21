<h1 align="center"><code>utop</code></h1>

<!-- badges: start -->
<p align="center">
  <a href="https://github.com/atsyplenkov/utop/actions/workflows/R-CMD-check.yaml"><img src="https://img.shields.io/github/actions/workflow/status/atsyplenkov/utop/R-CMD-check.yaml?style=flat&labelColor=1C2C2E&color=256bc0&logo=GitHub%20Actions&logoColor=white"></a>
  <a href="https://app.codecov.io/gh/atsyplenkov/utop"><img src="https://img.shields.io/codecov/c/gh/atsyplenkov/utop?style=flat&labelColor=1C2C2E&color=256bc0&logo=Codecov&logoColor=white"></a>
</p>
<!-- badges: end -->

> [!warning] 
> This package is still in development and API is subject to change.

`utop` is a nextgen fork and continuation of the original [`rtop` package](https://CRAN.R-project.org/package=rtop). It keeps the core top-kriging functionality for interpolation with variable spatial support while modernising the package internals and extending the feature set.

This fork is intended to be rewritten using more modern R tooling, including `mirai` for parallel processing, `S7` for a more explicit object system, and `stars` for spatiotemporal datasets in place of the deprecated `spacetime` package. It also aims to add broader support for universal kriging (hope so!).

## Installation

You can install the development version of `utop` from GitHub with:

``` r
# install.packages("pak")
pak::pak("atsyplenkov/utop")
```

## Testing

Run the regular test suite with:

``` r
devtools::test()
```

Some slower diagnostic integration tests are skipped by default. To include
them, set `UTOP_DIAGNOSTICS=true` before running tests, for example:

``` sh
UTOP_DIAGNOSTICS=true Rscript -e "devtools::test()"
```

## Linting and formatting

This project uses [`jarl`](https://github.com/etiennebacher/jarl) for linting and [`air`](https://github.com/posit-dev/air) for formatting R code. Configuration lives in `jarl.toml` and `air.toml` at the repository root.

```sh
jarl check R tests
air format R tests
```

## Roadmap

- [x] Migrate tests to `testthat`
- [x] Migrate to `roxygen2` documentation
- [ ] Replace spacetime with stars
- [ ] Add example data with timeseries
- [ ] Remove data.table and reshape2 deps
- [ ] Add mirai for parallel processing (?)
- [x] Add universal kriging support (!!)
- [ ] Add a proper OOP structure through S7 (?)

## Acknowledgements

`utop` builds on the original `rtop` package by Jon Olav Skøien and the method description published in Skøien et al. ([2014](https://doi.org/10.1016/j.cageo.2014.02.009)).

This project has made heavy use of AI-assisted pair programming (both `Kimi 2.6` and `Opus 4.7` via Opencode as of May 2026). It is highly doubtful that we would be able to put this together that quickly without AI help.

## License

`utop` is distributed under GPL-3, consistent with the `rtop`. For the upstream package lineage and original package reference, see the original [`rtop` CRAN page](https://CRAN.R-project.org/package=rtop).
