# multihostSIR

<!-- badges: start -->
[![R-CMD-check](https://github.com/kramera3/multihostSIR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/kramera3/multihostSIR/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE.md)
<!-- badges: end -->

An R package for building and simulating SIR (Susceptible-Infected-Recovered)
models with many interacting host species, where the per-species rates can be
set from body mass instead of being measured directly.

A lot of wildlife and zoonotic pathogens don't respect species boundaries.
SARS-CoV-2, avian influenza, rabies, and bovine TB all circulate in several
hosts at once. Standard one-host SIR models miss most of what makes those
systems interesting: transmission is asymmetric between species, contact is
shaped by who eats whom, and different hosts differ a lot in how susceptible
they are and how badly the infection hits them. This package is an attempt to
model those pieces without having to hand-write a new system of equations every
time you add a species.

Documentation comes in a few forms. There are two package vignettes (an
introduction and a fuller walkthrough of the model and study), which you can read
with `browseVignettes("multihostSIR")` once it's installed. The same material is
also available as plain markdown in [GUIDE.md](GUIDE.md) and as a runnable console
tour via `demo(getting_started, package = "multihostSIR")`, so nothing is lost if
your machine can't render the vignettes.

## Installing

The quickest way, which also checks your R version and pulls in anything
missing:

```r
source("https://raw.githubusercontent.com/kramera3/multihostSIR/main/install.R")
```

Or do it yourself:

```r
install.packages("remotes")
remotes::install_github("kramera3/multihostSIR")
```

Either `remotes` or `devtools` works. `remotes` is lighter, which is why it's
first here. By default neither builds the vignettes on install, which keeps the
install fast and avoids needing pandoc locally; if you want the vignettes built
and you have pandoc, add `build_vignettes = TRUE`:

```r
remotes::install_github("kramera3/multihostSIR", build_vignettes = TRUE)
```

The rendered vignettes are also on the package's GitHub page and on CRAN once
released, so you don't have to build them yourself to read them.

You need R 4.1 or later. The required packages (`deSolve`, `ggplot2`, `scales`,
`rlang`) install automatically. `openxlsx` is optional and only needed if you
want Excel output; `testthat` is only needed to run the tests.

Once it's installed, you can confirm everything is in place:

```r
library(multihostSIR)
check_dependencies()
#> multihostSIR dependency check
#>
#>   [ok  ] deSolve    >= 1.30    installed 1.40
#>   [ok  ] ggplot2    >= 3.4.0   installed 3.5.1
#>   [ok  ] scales     >= 1.2.0   installed 1.3.0
#>   [ok  ] rlang      >= 1.0.0   installed 1.1.4
#>   [MISS] openxlsx   any        not installed
#>   ...
```

If anything shows up as missing or out of date:

```r
install_dependencies()            # asks before it does anything
install_dependencies(ask = FALSE) # for scripts
```

## A first run

```r
library(multihostSIR)

set.seed(42)

# Draw a random community of 6 host species
comm <- build_community(6)
comm[, c("type", "mass", "N", "competence")]

# Build its transmission matrix and get the community-level R0
beta <- transmission_matrix(comm)
community_R0(comm, beta)

# Run a single epidemic and look at the per-species outcomes
sim <- simulate_epidemic(comm, duration = 365)
species_metrics(sim)

# Or the thing you'll actually use most: replicate many random communities
# across a range of species richness, timing each level as it goes
res <- run_multihost_sir(richness_levels = 3:8, n_iterations = 100, seed = 42)
summary(res)
plot(res, which = "r0")
```

If you'd rather read through a worked example than piece it together from the
snippet above, run the demo. It's the same material as GUIDE.md but as live code
you can step through:

```r
demo(getting_started, package = "multihostSIR")
```

## How users give input

Everything is driven through function arguments. There are no config files to
edit and nothing to set up globally. Three small objects hold the assumptions,
each with defaults, and you override whichever ones you care about:

```r
# 1. How body mass maps onto rates
allo <- allometry(
  a_recovery = 0.05    # slower recovery than the default 0.1
)

# 2. How random communities are drawn
spec <- community_spec(
  prey_prob               = 0.8,               # 80% of species are prey
  prey_masses             = c(0.02, 2, 50),    # body-mass pool to sample from
  predator_density_factor = 0.4                # predators less rare than default
)

# 3. Strength of each transmission route between trophic groups
coeffs <- trophic_coefficients(
  pred_to_prey = 0.1   # allow a route the default blocks
)

# Then hand them to the driver
res <- run_multihost_sir(
  richness_levels = 3:10,   # community sizes to simulate
  n_iterations    = 500,    # random communities per size
  duration        = 365,    # days per epidemic
  spec            = spec,
  allo            = allo,
  coefficients    = coeffs,
  seed            = 42,     # for reproducibility
  parallel        = TRUE,   # spread iterations across cores
  n_cores         = 4
)
```

If you just call `run_multihost_sir()` with no arguments beyond the richness
levels, you get the defaults, which are sensible starting points. Print any of
the three objects (`allo`, `spec`, `coeffs`) to see exactly what a run will do
before you start it.

The results come back as a list. `res$results` has one row per iteration with R0,
peak prevalence, timing to peak, community composition, and so on;
`res$trajectories` has the mean daily prevalence per richness level; `res$timing`
has the wall-clock time each level took. `summary(res)` collapses all of that
into a table by richness.

To save it out:

```r
export_results(res, dir = "output", formats = c("csv", "xlsx"))
```

## Bringing your own data

The run above generates communities at random and averages over many draws. If
you already know your system -- the actual species, their rates, and how they
transmit -- you can hand that in instead, and the model will use your numbers
rather than generating any.

A species table has one row per species. Give as much as you know:

| column | meaning | needed |
|---|---|---|
| `type` | "Prey" or "Predator" | always |
| `mass` | body mass, kg | always (unless every rate is given) |
| `mu` | mortality, per day | for a complete table |
| `gamma` | recovery, per day | for a complete table |
| `alpha` | disease mortality, per day | for a complete table |
| `N` | density, individuals per km2 | for a complete table |
| `competence` | within-species transmission | for a complete table |
| `species` | a name or label | optional |

Read it from a data frame, a CSV, or an Excel file:

```r
comm <- read_community("my_species.csv")
comm <- read_community("my_species.xlsx", sheet = "communities")
```

Any column you leave out is filled from body mass using the allometric
assumptions, exactly as the random generator would. Whatever you supply is kept
as given.

**When you supply everything, the answer is exact.** If your table is complete
*and* you also provide the transmission matrix, there is nothing random left to
average over, so the model runs once and returns the deterministic result -- no
iterations needed:

```r
comm <- read_community("my_species.csv")            # all rates supplied
beta <- read_transmission_matrix("my_matrix.csv")   # your own n-by-n matrix

fit <- run_fixed_community(comm, beta, duration = 365)
fit
fit$results        # one row: R0, peak prevalence, timing, ...
fit$species        # per-species outcomes
fit$trajectory     # the daily epidemic curve
```

If you supply a community but no matrix, the package builds the matrix from your
competence values and the trophic rules, and tells you it did so -- the run is
still a single pass, but it now rests on one built-in assumption. And if your
table is only partial, the missing pieces are generated, so you are back to
needing `run_multihost_sir()` with iterations to average over the randomness.

Every result prints an **"assumptions used"** block, so you can always see at a
glance which numbers were yours and which came from the built-in defaults.

Reading Excel needs the `readxl` package (`install.packages("readxl")`); CSV
needs nothing extra.

## The model

For species *i* in a community of *n*:

```
dS_i/dt = b_i - mu_i * S_i - Lambda_i * S_i
dI_i/dt = Lambda_i * S_i - (mu_i + gamma_i + alpha_i) * I_i
dR_i/dt = gamma_i * I_i - mu_i * R_i
dC_i/dt = Lambda_i * S_i
```

`Lambda_i = sum_j beta[i,j] * I_j` is the force of infection on species *i*, and
`C_i` just accumulates cumulative incidence (it doesn't feed back into anything).
Births are held constant at `b_i = mu_i * N_i`, so the disease-free equilibrium
sits at `S_i = N_i`. One consequence worth flagging: total population size is
*not* conserved once virulence `alpha_i` is positive, which is intentional:
disease-induced mortality is allowed to draw down abundance.

### Where the rates come from

If you don't have measured rates for a species, they're derived from its body
mass `M` (in kg) via `Y = a * M^b`. The time unit is one day throughout.

| Rate | Meaning | Default | Why this exponent |
|---|---|---|---|
| mu | background mortality (per day) | `0.005 * M^-0.25` | Physiological rates and lifespan scale as roughly `M^-0.25`; bigger animals live longer (Kleiber 1932; Peters 1983). |
| gamma | recovery rate (per day) | `0.1 * M^-0.25` | Process durations go as `M^0.25`, so rates go as `M^-0.25`; bigger animals recover slower (Lindstedt & Calder 1981). |
| N | carrying capacity (ind. per km2) | `500 * M^-0.75` | Damuth's law; bigger animals are rarer (Damuth 1981). |
| alpha | disease-induced mortality (per day) | `0.01 * M^-0.25` | Same quarter-power scaling. |

### Transmission between species

The diagonal of the beta matrix is the within-species rate,
`beta_ii = competence_i / N_i`. Off-diagonal entries are the geometric mean of
the two within-species rates, `sqrt(beta_ii * beta_jj)`, scaled by a coefficient
that depends on the trophic types involved. Rows are the recipient of infection,
columns are the source. The defaults are:

| recipient \ source | Prey | Predator |
|---|---|---|
| **Prey** | 0.50 | 0.00 |
| **Predator** | 0.25 | 0.25 |

So under the defaults predators pick up infection from their prey but never pass
it back, which shows up as a block of zeros in the matrix. Change any of these
with `trophic_coefficients()`.

The community R0 is the dominant eigenvalue of the next-generation matrix
`K[i,j] = beta[i,j] * N[i] / (gamma_j + mu_j + alpha_j)`.

## What's in the package

The functions group roughly like this:

- Setup: `check_dependencies()`, `install_dependencies()`
- Defining assumptions: `allometry()`, `community_spec()`,
  `trophic_coefficients()`
- Building a community: `build_community()`, `transmission_matrix()`
- Bringing your own data: `read_community()`, `read_transmission_matrix()`,
  `run_fixed_community()`
- R0: `next_generation_matrix()`, `community_R0()`
- Simulating: `simulate_epidemic()`, `species_metrics()`, `simulate_once()`
- The main driver: `run_multihost_sir()`, with `summary()`,
  `timing_complexity()`, `predator_summary()`, `competence_summary()`
- Plots, by richness: `plot_r0_richness()`, `plot_peak_prevalence()`,
  `plot_cumulative_richness()`, `plot_days_to_peak_richness()`,
  `plot_trajectories()`
- Plots, macro-drivers: `plot_r0_trophic()`, `plot_r0_body_mass()`,
  `plot_r0_density()`, `plot_r0_competence()`, `plot_cumulative_competence()`,
  `plot_peak_density()`, `plot_days_to_peak_r0()`
- Plots, predator structure: `plot_peak_predator()`, `plot_r0_predator()`,
  `plot_trophic_partition()`
- Plots, methods: `plot_transmission_matrix()`, `plot_allometry()`,
  `plot_timing()` (add `log_scale = TRUE` for the log-log timing)
- All of the above at once: `plot_gallery()` returns them in a named list
- Fixed-run plots: `plot_fixed_trajectory()`, and `plot()` on a fixed run
- Saving: `export_results()`

Every exported function has help with runnable examples, so `?run_multihost_sir`
and friends are the fastest reference.

## Repository layout

```
multihostSIR/
  install.R                bootstrap installer (checks deps, installs, verifies)
  GUIDE.md                 the full walkthrough (plain-markdown mirror)
  R/                       package source, split by topic
  man/                     help pages, generated by roxygen2 (don't edit by hand)
  vignettes/               the two package vignettes (.Rmd)
  demo/getting_started.R   the guided tour as runnable code
  tests/testthat/          unit tests
  analysis/                the full simulation study (not part of the package)
  DESCRIPTION, NAMESPACE    package metadata
```

`analysis/`, `install.R`, and `GUIDE.md` are in `.Rbuildignore`, so they live in
the repo but aren't part of the installed package.

The full study is a plain R script, not an Rmd:

```bash
Rscript analysis/reproduce_study.R
```

It writes figures to `analysis/figures/`, tables to `analysis/output/`, and a
console log alongside them.

## Reproducibility

Same seed, same results:

```r
a <- run_multihost_sir(3:5, n_iterations = 50, seed = 42)
b <- run_multihost_sir(3:5, n_iterations = 50, seed = 42)
identical(a$results$R0, b$results$R0)   # TRUE
```

One thing to know if you go parallel: a serial run and a parallel run with the
same seed won't match, because the parallel path uses L'Ecuyer streams
(`clusterSetRNGStream`) rather than the default RNG. Within a mode it's
deterministic. So if you need bit-for-bit identical output, keep `parallel`,
`n_cores`, and `seed` all fixed. The seed is set after an untimed warm-up batch,
so warm-up never shifts your results.

## Developing

```r
install.packages(c("roxygen2", "testthat", "knitr", "rmarkdown"))

roxygen2::roxygenise()    # regenerate man/ and NAMESPACE from the #' comments
testthat::test_local()    # run the tests
```

If you edit any of the `#'` documentation blocks, run `roxygenise()` again before
committing so `man/` stays in sync. To preview a vignette while editing it,
`rmarkdown::render("vignettes/multihostSIR.Rmd")` (needs pandoc). See
[CONTRIBUTING.md](CONTRIBUTING.md) if you want to send a change.

## Citing

```r
citation("multihostSIR")
```

## References

Damuth, J. (1981) Population density and body size in mammals. *Nature* 290, 699-700.

Diekmann, O., Heesterbeek, J.A.P. & Metz, J.A.J. (1990) On the definition and the computation of the basic reproduction ratio R0 in models for infectious diseases in heterogeneous populations. *Journal of Mathematical Biology* 28, 365-382.

Kleiber, M. (1932) Body size and metabolism. *Hilgardia* 6, 315-353.

Lindstedt, S.L. & Calder, W.A. (1981) Body size, physiological time, and longevity of homeothermic animals. *Quarterly Review of Biology* 56, 1-16.

Peters, R.H. (1983) *The Ecological Implications of Body Size*. Cambridge University Press.

## License

MIT. See [LICENSE.md](LICENSE.md).
