_default:
    just --list

# Format R code
fmt:
    air format .

# Run jarl linter on the R package
lint-r:
    jarl check .

# Alias for lint-r
lint: lint-r

# Regenerate Rd files and NAMESPACE from roxygen comments
docs:
    Rscript -e "roxygen2::roxygenise()"

# Run testthat tests interactively
test:
    Rscript -e "devtools::test()"

# Build and check the R package as CRAN would
check:
    #!/usr/bin/env bash
    set -euxo pipefail
    rm -rf utop_*.tar.gz utop.Rcheck
    R CMD build .
    R CMD check --as-cran utop_*.tar.gz
    rm -rf utop_*.tar.gz utop.Rcheck

# Build and check the R package without building the manual
check-no-manual:
    #!/usr/bin/env bash
    set -euxo pipefail
    rm -rf utop_*.tar.gz utop.Rcheck
    R CMD build .
    R CMD check --no-manual utop_*.tar.gz
    rm -rf utop_*.tar.gz utop.Rcheck

# Build and install the R package locally
install:
    #!/usr/bin/env bash
    set -euxo pipefail
    rm -rf utop_*.tar.gz
    R CMD build .
    R CMD INSTALL utop_*.tar.gz

# Remove build artifacts
clean:
    rm -rf utop_*.tar.gz utop.Rcheck
