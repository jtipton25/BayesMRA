# BayesMRA

An R package for Bayesian spatial and spatio-temporal regression using sparse multi-resolution matrices. This package includes functions for setting up the multi-resolution process and a function `mcmc_mra()` for fitting a spatial model using Markov Chain Monte Carlo.


## Installation

## Package Installation

If this is your first time installing this package, [make sure you have the appropriate compiler](#make-sure-you-have-the-appropriate-compiler).


### Installation from gitHub using `devtools`

To install the package from gitHub, [make sure you have the appropriate compiler](#make-sure-you-have-the-appropriate-compiler). Next, make sure you have the `devtools` packge installed. You can install the `devtools` package using 

```r
install.packages("devtools")
```

Once the `devtools` library is installed, you can install the `BayesMRA` library using

```r
devtools::install_github("jtipton25/BayesMRA")
```
**Note: It is recommended to reguarly check for updates by using `devtools::install_github("jtipton25/BayesMRA")` frequently as this package is in active development and regularly undergoes large changes**

### Installation from CRAN

Installation from CRAN is not currently supported (hopefully coming soon). 

<!-- 
Therefore
```r
install.packages("BayesMRA")
``` 
will fail.
-->


### Make sure you have the appropriate compiler

- **Linux**: It is assumed on Linux systems, the appropriate compiler is available

- **Windows**: For Windows, download and install RTools. For R (>=4.0.0) follow the instructions [here](https://cran.r-project.org/bin/windows/Rtools/) -- older versions of R follow the instructions [here](https://cran.r-project.org/bin/windows/Rtools/history.html).

- **MacOS**: If you are on a Mac, make sure you have an openMP supported compiler -- see [here](https://thecoatlessprofessor.com/programming/cpp/r-compiler-tools-for-rcpp-on-macos/) for instructions on how to get this setup. Follow the instructions for your specific version of R


