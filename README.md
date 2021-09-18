# Metrics Paper Figures #

Code and data for generating data figures used in the main and supplemental sections of the manuscript *"A framework for evaluating edited cell libraries created by massively parallel genome engineering"*

## Overview ##

* [Quickstart](#markdown-header-quickstart)
* [Organization of contents](#markdown-header-organization)

## Quickstart ##

* Create the `cell_lib_eval_paper` conda environment: `conda env create -f environment.yaml`
* Activate conda environment: `conda activate cell_lib_eval_paper`
* Generate the figures used in the main section of the paper: `R -e "rmarkdown::render('Rmd/main_figures.Rmd', output_file='../main_figures.html')"`
* Generate the figures used in the supplemental section of the paper: `R -e "rmarkdown::render('Rmd/supp_figures.Rmd', output_file='../supp_figures.html')"`

## Organization ##

* Inputs
    * `README.md` - this file
    * [`environment.yaml`](environment.yaml) - definition of `cell_lib_eval_paper` conda environment required to generate figures
    * [`Rmd/`](Rmd/) - Rmarkdown code used to generate figures
    * [`R/`](R/) - common R code used by the Rmarkdown files in `Rmd/`
    * [`Data/`](Data/) - data used in the generation of figures
* Outputs
    * [`png/`](png/) - directory with individual png files of figures and figure sub-panels
    * [`main_figures.html`](main_figures.html) - figures used in the main section of the paper
    * [`supp_figures.html`](supp_figures.html) - figures used in the supplemental section of the paper

## Deeper Dive ##

Figure 4 from the supplemental section is generated based on simulated data.  For convenience, the simulated data has been pre-generated and stored in `Data/rdata/richness_sim.RData`.  If desired, reproduction of the simulated data is possible by following these additional steps.

* Install additional required packages, either via conda as described below, or from CRAN or with `install.packages()` in an R session
    * `conda install -c conda-forge r-doparallel`
    * `conda install -c conda-forge r-foreach`
* Run the simulation to create results in `Data/rdata/richness_sim.RData`.  This will benefit from being run on a multi-core machine - on a system with 56 cores it took just under two hours.
    * `R CMD BATCH --vanilla R/generate_sim.R`
