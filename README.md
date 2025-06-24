# rmi
Tools for analyzing and plotting RMI data.

Author: Lori Garzio (lgarzio@marine.rutgers.edu)

Funding provided by New Jersey’s Research & Monitoring Initiative (RMI) (New Jersey Department of Environmental Protection, New Jersey Board of Public Utilities)

## Installation Instructions
Add the channel conda-forge to your .condarc. You can find out more about conda-forge from their website: https://conda-forge.org/

`conda config --add channels conda-forge`

Clone the rmi repository

`git clone https://github.com/lgarzio/rmi.git`

Change your current working directory to the location that you downloaded rmi. 

`cd /Users/garzio/Documents/repo/rmi/`

Create conda environment from the included environment.yml file:

`conda env create -f environment.yml`

Once the environment is done building, activate the environment:

`conda activate rmi`

Install the toolbox to the conda environment from the root directory of the rmi toolbox:

`pip install .`

The toolbox should now be installed to your conda environment.

## Steps
1. Using the rucool (dataset_archiving repo)[https://github.com/rucool/dataset_archiving/tree/master], download and QC glider datasets.
