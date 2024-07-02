## Reproduction Package

The scripts in this folder are used for generating
the plots and processing the data used in the paper.
We will describe the different steps of the data processing
pipeline and provide bash scripts to reproduce the pipeline.

Whenever you use scripts from this reproduction package, please cite this Zenodo artifact and the published article.

There are two different ways of using this package:
1) Directly via Python on your local 
machine.
2) In a docker container (**recommended, if you are not familiar with executing 
scripts**).

### WARNING:
Some of the data processing conducted on the PetDB data set requires advanced computational 
power and otherwise will be extremely time-expensive. Therefore, we also added pre-computed data 
files that will shorten the computation time. You 
can either use those pre-computed files or run all computations on your own.

<span style="color: #FF0000">If you decide to
run all computations on your own, this will over-write the pre-computed data!</span>

If you want to use the pre computed data, use `./run_pipeline_pre_computed.sh`at the respective step.

If you want to reduce computation even further, you can exclude the Monte Carlo Simulations for the PetDB data and use `./run_pipeline_noMCS.sh`

###  Option 1: Setting up the environment
The scripts shown here requires a bunch of python packages to work well.
Among them is pyrolite, pandas, numpy, etc..
We recommend to setup a new environement with [miniconda](https://docs.anaconda.com/free/miniconda/index.html).
 `conda env create -f environment.yml` should do the job on your machine and create 
an environment _paper_ for reproducing these scripts. 

### Simulating NAA and ID Data
We start with a data from the PetDB database as explained in the paper.
From this data, we create a an NAA and a ID data set using the `generate_test_data.py` script.

### Computing the Models.
The `lambdatau_fitting.py` script and the `polyfit_fitting.py` scripts
are used to compute models for the original data.

### Reproduce Figures

To reproduce figures, fun the `./run_pipeline.sh` command.
Otherwise, run commands for each figure on its own.
The `run_pipeline.sh` descripes how to invoke the script for each figure. 


### Option 2: Running in the docker container
#### Build the container

Just run `docker build -t paper_ree --platform=linux/amd64 .`.


#### Use the container

Run `docker run -v ${PWD}/figures:/figures -it paper_ree`.
Inside the container run:
- `conda activate paper`
- `./run_pipeline.sh`
- find results in the output folder 