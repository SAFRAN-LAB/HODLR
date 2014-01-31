# HODLR Examples

## Kepler

This is a test on real data from the NASA Kepler mission. The test is
implemented in `./kepler_test.cpp` and the data are given in
`./kepler-10-sc.csv`. The data are "short cadence" observations of the star
Kepler 10 and the columns are time (measured in KBJD), relative flux (in
arbitrary units), and the uncertainties on the fluxes. To run this example,
after `make`-ing while still in the `build` directory, run
```
bin/kepler < ../examples/kepler-10-sc.csv
```
