# treeme.R

R/ggtree script for automated tree plotting

### Install

mamba create -n treeme r-base r-optparse r-ggplot2 r-cowplot r-dplyr r-pacman bioconductor-treeio bioconductor-ggtree

### Testing

Testing uses `bats-core`. To run tests:

```
bats test/test.bat
```