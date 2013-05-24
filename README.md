## Usage

To replicate the differential expression analysis shown in Filone 2013, clone (or [download](https://github.com/nachocab/vaccinia_filone_2013/archive/master.zip) and extract) the current repository.

Then, launch an R session and run the following:
```
# set the current working directory to the location of the current repository in your local drive
setwd("path/to/vaccinia_filone_2013-master")

# install required R packages
install.packages(c("knitr", "markdown", "edgeR", "stringr", "bvenn", "RColorBrewer", "ggplot2", "MASS", "scales", "grid", "reshape2"))

# load required packages
library(knitr)
library(markdown)

knit2html("analysis.Rmd")
browseURL("analysis.html")
```
