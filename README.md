# confoundr
This preliminary R package implements the methods described in <a href="https://www.ncbi.nlm.nih.gov/pubmed/?term=jackson+jw+diagnostics">Diagnostics for Confounding of Time-Varying and Other Joint Exposures. Epidemiology 2016;27(6):859-69"</a>.

To install the package, use the following code:

```
install.packages("devtools",dependencies=TRUE)
library(devtools)
install_github("jwjackson/confoundr",dependencies=TRUE,build_vignettes=TRUE)
```

To load the package, use the following code:

```
library(confoundr)
```

The package now contains the documentation and a vignette. However, if you wish to use the functions directly, download the file Rfunctions_1_0_12.r' in the R directory. A PDF manual can be found in the INST directory. The toy data 'example_sml.csv' used in the manual's example can be found in the R directory.

For questions please contact me

John W. Jackson, ScD
