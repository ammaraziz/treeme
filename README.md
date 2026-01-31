# treeme.R

R/ggtree script for automated tree plotting

### Install

```
mamba create -n treeme r-base r-optparse r-ggplot2 r-cowplot r-dplyr r-ggnewscale r-rlang bioconductor-treeio bioconductor-ggtree=4.0.4
```

### Run
```
Rscript treeme.R
```

### Testing

Testing uses `bats-core`. To run tests:

```
bats test/test.bat
```