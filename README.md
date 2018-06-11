This repository contains the data and source code for the following paper:

* J. Urbano and T. Nagler, "[Stochastic Simulation of Test Collections: Evaluation Scores](http://julian-urbano.info/files/publications/065-stochastic-simulation-test-collections-evaluation-scores.pdf)", *International ACM SIGIR Conference on Research and Development in Information Retrieval*, 2018.

A [single ZIP file](https://github.com/julian-urbano/sigir2018-simulation/archive/master.zip) can be downloaded as well.

## Project Structure

* `data/` Input data files.
* `output/` Generated output files.
* `src/` Source code in R.
* `scratch/` Temporary files generated in the process.

All code is written for [R](https://www.r-project.org). You will need the following packages installed from CRAN: `simIReff`, `VineCopula`, `rvinecopulib`, `doParallel`.

## How to reproduce the results in the paper 

The source files in `src/` need to be run in order. You can run each file individually by running `Rscript src/<file>.R`. They will store intermediate data in `scratch/` and the final data in `output/`.

**It is important that you always run from the base directory**.

1. `src/01-margins.R` fits all the marginal distributions (section 4.1).
2. `src/02-transform.R` transforms all distributtions fitted in step 1 to have the same expected value as in the original data (section 4.1).
3. `src/03-margins_select.R` selects the best distributions from step 1 (section 4.1).
4. `src/04-transform_select.R` selects the best distributions from step 2 (section 4.2).
5. `src/05-bicop_types.R` fits all possible bivariate copulas (section 4.2).
6. `src/06-cop_fit.R` fits all the full gaussian and R-vine copulas (section 4.2).
7. `src/07-cop_sim.R` simulates from the full copulas in step 6 and records observed means and variances (section 4.3).
8. `src/10-power.R` runs the first sample application (section 5.1).
9. `src/11-type_I.R` runs the second sample application (section 5.2).
7. `src/99-paper.sh` generates all figures in the paper and stores them in `output/figs/`.

It takes several days to run all the code, so it is ready to run in parallel. Most of the above code parallelizes using function `foreach` in R's package [`doParallel`](https://cran.r-project.org/web/packages/doParallel/index.html). In particular, it will use all available cores in the machine except one. Check the individual code files to modify this behavior.

## Custom test collections

You can easily run the code with your own initial test collection. Add the matrix of topic-by-system scores in `data/` using the name `<collection>_<measure>.csv` (see for instance file [`data/web2012_ap.csv`](/data/web2012_ap.csv)). Then, edit file `src/common.R` to add the new data:

```r
.COLLECTIONS <- c("web2010", "web2011", "web2012", "web2013", "web2014")
.MEASURES <- c("ap", "ndcg20", "err20", "p10", "p20", "rr")
```

Note that the code will run for all combinations of collection and measure. For more specific modifications, edit the corresponding source file in `src/` (see above). Note also that the script `src/99-paper.R` is only intended to generate the figures in the paper. If you customize something and want a similar analysis, you will need to extend this script yourself.

## How to simulate your own data

A full R package has been developed for this: [`simIReff`](https://cran.r-project.org/web/packages/simIReff/index.html). Please refer to it's documentation for details. The simplest use case would be something like:

```r
d <- read.csv("data/web2012_p10.csv") # read original data
effs <- effDiscFitAndSelect(d, support("p10")) # fit and select margins
cop <- effcopFit(d, effs) # fit copula

y <- reffcop(1000, cop) # simulate new 1000 topics
# compare observed vs. expected mean
expected <- sapply(effs, function(e) e$mean)
observed <- colMeans(y)
plot(expected, observed)
```

## License

* The TREC results in `data/` are anonymized and posted here with permission from the organizers.
* Databases and their contents are distributed under the terms of the [Creative Commons Attribution-ShareAlike 4.0 International License](http://creativecommons.org/licenses/by-sa/4.0/).
* Software is distributed under the terms of the [MIT License](https://opensource.org/licenses/MIT).

When using this archive, please [cite](CITE.bib) the above paper:

    @inproceedings{urbano2018simulation,
      author = {Urbano, Juli\'{a}n and Nagler, Thomas},
      booktitle = {International ACM SIGIR Conference on Research and Development in Information Retrieval},
      title = {{Stochastic Simulation of Test Collections: Evaluation Scores}},
      year = {2018}
    }
