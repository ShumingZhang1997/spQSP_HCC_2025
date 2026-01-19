# Spatial QSP (HCC) Model

## Ubuntu Operating System Configuration (Only required for Windows User)

This step helps to setup the Ubuntu operating system in for Windows user via virtual machine.
1. The virtual machine host: VirtualBox is available at https://www.virtualbox.org/
2 . The Ubuntu Desktop image (Latest version 18.04) is available at   [https://releases.ubuntu.com/18.04/]
3. Enter the “Oracle VM VirtualBox Manager”, press “New” to create the virtual machine with all default settings. (Recommend allocate 20 GB for storage and 2 GB of RAM)

**Notice**: All following operations should be done in the Linux operating systems (the virtual machine), NOT Windows.

## Required Library Installation
Libraries: **SUNDIALS**: version:4.0.1; **Boost**: version 1.70.0

-SUNDIALS <br />
1. Download is available at: https://computing.llnl.gov/projects/sundials/sundials-software <br />
The following files are downloaded: 
 `sundials-4.0.1.tar.gz`

2. Decompress Archieve: <br />
`$ tar xzf sundials-4.0.1.tar.gz`

3. Install cmake if not already available: <br />
`$ sudo apt install cmake-curses-gui`

4.	Create install and build directories: <br />
`$ mkdir -p ~/lib/sundials-4.0.1` <br />
`$ mkdir -p ~/Downloads/sundials-build` <br />
`$ cd ~/Downloads/sundials-build` <br />

5. Configuration <br />
`$ ccmake ~/Downloads/sundials-4.0.1` <br />
Press c key to enter configuration interface
Set install directory: CMAKE_INSTALL_PREFIX set to `~/lib/sundials-4.0.1`
Set example install directory: EXAMPLE_INSTALL_PATH set to `~/lib/sundials-4.0.1/examples`
Press c repeatedly to process configuration; press g to generate Makefile and exit configuration.
6. Build <br />
From `~/Downloads/sundials-build/` <br />
`$ make` <br />
`$ make install`

-Boost Version 1.70.0 <br />

1. Source code available at: https://www.boost.org/users/history/version_1_70_0.html <br />
The following files are downloaded: <br />
`boost_1_70_0.tar.gz`

2.	Decompress the archive: <br />
`$ tar xzf boost_1_70_0.tar.gz` <br />

Official instructions is available at:
https://www.boost.org/doc/libs/1_70_0/more/getting_started/unix-variants.html

3. Building separately-compiled boost libraries <br />
`$ cd ~/Downloads/boost_1_70_0` <br />
`$ ./bootstrap.sh --prefix=$HOME/lib/boost_1_70_0` <br />
`$ ./b2 install` <br />

## Model Simulation 
The Makefile of this model is available at: `~/HCC/HCC_single/linux/` <br />

To prepare spQSP for a simulation, write:
`$ make HCC_s_sim` <br />

Then, to show all options to configure the simulation:
`$ ./HCC_s_sim -h` <br />

### Command-line options

The executable supports the following options:

| Option | Default | Description |
|--------|---------|-------------|
| `-h, --help` | — | Produce help message |
| `-s, --seed` | `0` | Seed value |
| `-t, --time` | `0` | Total number of simulation steps |
| `-p, --param-file` | — | Parameter file name |
| `-o, --output-path` | `defaultOut` | Output file base path |
| `--outParam` | `outParam.xml` | Save a copy of parameter file |
| `-B, --brief` | — | Print brief tracking info to stdout |
| `-S, --stats` | — | Print simulation statistics |
| `--stats-interval` | `1` | Interval to save stats |
| `-G, --grid` | `0` | Grid print mode (0: none, 1: cell, 2: grid, 3: both) |
| `--grid-interval` | `1` | Interval to print grid information |
| `--save-state-start` | — | Time step to start saving state |
| `--save-state-interval` | — | Interval for saving state |
| `--load-state` | — | Load saved state file |

---

### Example command

Run a simulation and save outputs to `Outputs/`:

`$ ./HCC_s_sim -t 280 -p ../resource/param_all_test.xml -o Outputs -B -S -G 1`
It is recommended to run `$ clean` <br /> before running `$ make HCC_s_sim` <br /> if the code has been modified in between simulations.

# Model Calibration Framework

This repository implements a model calibration workflow using Approximate Bayesian Computation Sequential Monte Carlo (ABC-SMC) with **pyABC**.

## Environment

The calibration framework has been tested with the following environment:

- **Python**: 3.9.18  
- **pyABC**: 0.12.13  
- **Dask**: 2024.7.1  
- **dask-jobqueue**: 0.8.5
- **R**: 4.3.0 *(required for integrated R scripts and analyses)*

These versions ensure compatibility for distributed and parallel computation during model calibration.

## Overview

This pipeline enables:

- Bayesian parameter inference via ABC-SMC  
- Parallelized model evaluations using Dask  
- Scalable HPC execution with `dask-jobqueue`

## Run Calibration
To run the calibration framework in the paper:
`$ python model_calibration_spqsp.py` <br />

## Citation
> Schälte, Y., Klinger, E., Alamoudi, E., & Hasenauer, J. (2022). *pyABC: Efficient and robust easy-to-use approximate Bayesian computation*. **Journal of Open Source Software, 7**(74), 4304. https://doi.org/10.21105/joss.04304
