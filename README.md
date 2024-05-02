# Last update: 02/05/2024
## Benedetti et al. (in prep. for Global Change Biology) - Global copepod functional diversity patterns

This repository contains R scripts developed during my postdoc at ETH Zürich, D-USYS, IBP, UP group.
The present R scripts were developed to produce the results of our study focusing on the global contemporary and future patterns of marine copepod functional diversity, and how it relates to the functioning of the biological carbon pump.

- Scripts labelled as 'RSCRIPTBATCH' are those that were run on a local cluster to perform those functions defined within the script in parallel using the 'parallel' R package (R Core Team, 2020) within a mclapply().
- Scripts that were given a number correspond to those where I formatted and mined the main datasets involved. They usually contain organized sequences of code where data are being read, examined, reformatted, analyzed and plotted. The number given to a script mainly corresponds to a temporal marker, and not necessarily a direct sequence from one numbered script to another.
- The end of the scripts' names usually indicate their main purposes. A more detailed list of the scripts' goals and content is usually given in the beginning of the numbered scripts.
