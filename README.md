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

- **`/Source_codes/`**  
  Contains the source code of the modelers used in the experiments.

- **`/Experiments/`**  
  Contains the `SConstruct` file responsible for generating the modeled data used in the comparisons.

- **`/scripts_py/`**  
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

### 2. Install Madagascar

    Please refer to the [official MADAGASCAR package documentation](https://ahay.org/wiki/Main_Page) for instructions on [downloading](https://ahay.org/wiki/Download) and [installing](https://ahay.org/wiki/Installation) the package.
