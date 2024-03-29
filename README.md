<!-- badges: start -->
[![R-CMD-check](https://github.com/jepusto/ARPobservation/workflows/R-CMD-check/badge.svg)](https://github.com/jepusto/ARPobservation/actions)
[![Codecov Status](https://codecov.io/gh/jepusto/ARPobservation/branch/master/graph/badge.svg)](https://app.codecov.io/github/jepusto/ARPobservation?branch=master)
[![](http://www.r-pkg.org/badges/version/ARPobservation)](https://CRAN.R-project.org/package=ARPobservation)
[![](http://cranlogs.r-pkg.org/badges/grand-total/ARPobservation)](https://CRAN.R-project.org/package=ARPobservation)
[![](http://cranlogs.r-pkg.org/badges/last-month/ARPobservation)](https://CRAN.R-project.org/package=ARPobservation)
<!-- badges: end -->

ARPobservation provides a set of tools for simulating the sort of data generated during direct observation of a behavior over time. It works by first simulating a behavior stream based on an alternating renewal process, given specified distributions of event durations and interim times. The following distributions are currently supported:

* Exponential
* Gamma with fixed shape parameter
* Mixture of two gammas
* Weibull with fixed shape parameter
* Uniform(0,x)
* Constant (degenerate)

Different methods of recording data can then be applied to the simulated behavior stream. Functions are provided for the following recording methods: 

* continuous duration recording 
* event counting
* momentary time sampling
* partial interval recording
* whole interval recording

