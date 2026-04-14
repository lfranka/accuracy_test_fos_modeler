# Accuracy Validation and Numerical Dispersion Analysis of a 2D Vector Acoustic Finite-Difference Staggered-Grid Time-Domain Modeler

## Overview
This repository contains the source code and experimental scripts associated with the paper:

> *Accuracy Validation and Numerical Dispersion Analysis of a 2D Vector Acoustic Finite-Difference Staggered-Grid Time-Domain Modeler*  
> Lucas K. F. Queiroz and Daniel Leal Macedo

The main goal of this work is to establish a **pragmatic workflow for validating the accuracy of a numerical modeling tool** through three numerical experiments:

1. **Modeler vs. Semi-Analytical Solution**  
   Comparison between modeled data and semi-analytical solutions under variations of modeling parameters (e.g., grid spacing and frequency).

2. **Group Velocity Measurement**  
   Comparison between theoretical group velocity (derived from the discretized equations) and numerically measured values.

3. **Phase Velocity Measurement**  
   Comparison between theoretical phase velocity (derived from the discretized equations) and numerically measured values.

---

## Repository Structure

- **`/source_code/`**  
  Contains the source code of the modelers used in the experiments.

- **`/experiments/`**  
  Contains the `SConstruct` file responsible for generating the modeled data used in the comparisons.

- **`/experiments/scripts_py/`**  
  Python scripts implementing the three experiments:
  - `comp_data.py` → Modeler vs. semi-analytical comparison  
  - `med_cg.py` → Group velocity measurement  
  - `med_cp.py` → Phase velocity measurement  

---

## Installation Guide

### 1. Clone the Repository
```bash
git clone https://github.com/lfranka/accuracy_test_fos_modeler.git
cd accuracy_test_fos_modeler
```
### 2. Install Madagascar

    Please refer to the [official MADAGASCAR package documentation](https://ahay.org/wiki/Main_Page) for instructions on [downloading](https://ahay.org/wiki/Download) and [installing](https://ahay.org/wiki/Installation) the package.

### 3. Install the Source Code

   Once Madagascar is installed, navigate to the users folder and create a new directory with your name. 
   ```bash
   cd madagascar-4.2/user
   mkdir your_name 
   ```
   Then copy all the files from the sources codes folder to your directory in the users directory.
   ```bash
   cp -v -r sources_code/* madagascar-4.2/user/your_name  
   ```
   After that run scons in your user folder
   ```bash
   cd madagascar-4.2/user/your_name
   scons
   ```
   And finally run scons in the madagascar folder
   ```bash
   cd ../..
   sudo scons install
   ```

## Runing the Experiments

   Go to the experiments directory and run scons
   ```bash
   cd experiments
   scons
   ```
   Then go to scripts_py directory and run the three scripts to see and save the results
   ```bash
   python name_scripts.py
   ```
   The data for group and phase velocity is save in data directory and the figures of the results is save in figs_comp figs directory.

## Abstract
Seismic numerical modeling plays a crucial role in subsurface imaging and resource exploration, with finite difference methods widely employed to simulate wave propagation in complex media. In this study we provide tools for evaluating the accuracy of a 2D vector-acousticc modeler. To this end we conducted three experiments: (i) comparison between modeled and semi-analytical seismograms to assess waveform fidelity; (ii) cross-correlation analysis to estimate group velocity and compare it with theoretical values; and (iii) Fourier analysis of two modeled data to estimate phase velocity and, again, compare it with the theoretical one. These tests aim to quantify the modeler's accuracy and stability, providing insights into its suitability for seismic modeling applications.

## Citation
If you use this code or dataset in your research, please cite the paper:

**Not published yet**
