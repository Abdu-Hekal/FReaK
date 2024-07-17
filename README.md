# FReaK - Falsification using Reachability and Koopman 

## Overview

A MATLAB toolkit for falsification using Koopman surrogate models. 

## Table of Contents

- [Introduction](#introduction)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Usage](#usage)
- [Citing FReaK](#Citing-FReaK)
- [Contact](#contact)

## Introduction

FReaK is a novel framework that generates counterexamples for falsification of black-box models.
We employ an iterative refinement technique based on surrogate modelling. In particular, simulations are used 
to construct a surrogate model for the system dynamics using data-driven Koopman operator linearization. 
The reachable set of states are then computed and combined with an encoding of the signal temporal logic specification in a mixed-integer linear program (MILP).
To determine the next sample, the MILP solver computes the least robust trajectory inside the reachable set of the surrogate model.
The trajectory's initial state and input signal are then executed on the original black-box system, where the specification is either falsified or additional simulation data is generated that we use to retrain the surrogate Koopman model and repeat the process.

## Prerequisites

FReaK requires the following prerequisites:

- **MATLAB** - You can download MATLAB from [MathWorks website](https://www.mathworks.com/).
- **Breach** - Download Breach toolbox from github repository [here](https://github.com/decyphir/breach/tree/master)
- **CORA(v2022)** - Download CORA v2022 from the archive page on the official [CORA website](https://tumcps.github.io/CORA/). Follow the instructions in the corresponding manual for installation.
- **YALMIP** - See CORA2022 manual for best installation description
- **Python(>=3.8)** - Install Python version 3.8 or higher.
- **AutoKoopman** - Install AutoKoopman directly from PyPI using: `pip install autokoopman`
- **Gurobi Optimization Toolbox** - We use Gurobi solver for backend optimization, though other solvers may be used. Follow the instructions provided on the [Gurobi website](https://support.gurobi.com/hc/en-us/articles/4533938303505-How-do-I-install-Gurobi-for-Matlab-) to download and install Gurobi for matlab

## Installation

We provide a setup script to help with installation. The script automatically installs prerequisites where possible 
or outlines which are missing via the command window. The script also makes necassary modifications to external toolboxes
and removes any conflicting files from path. Add the FReaK folder and subfolders to path and run setup file
**setupKF.m**

```matlab
% Run setup file
setupKF
```

The script also ensure that only necessary prerequisites for our framework for external toolboxes (like CORA) are installed.
For full functionality of each external toolbox, refer to their corresponding instructions.

Note that the only prerequisites that need manual installation are:
- **MATLAB**
- **Python(>=3.8)**
- **Symbolic Math Toolbox** - A MATLAB toolbox required by CORA. Install from MATLAB Add-On explorer, see [here](https://uk.mathworks.com/help/matlab/matlab_env/get-add-ons.html)
- **Gurobi Optimization Toolbox** - Download from the Gurobi website [Gurobi website](https://support.gurobi.com/hc/en-us/articles/4533938303505-How-do-I-install-Gurobi-for-Matlab-)
and requires a license

## Usage

To test successful installation, run example script that falsifies a Vanderpol model and plots the falsification trace:

```matlab
% Run Vanderpol falsification 
falsifyVanderpol
```



## Citing FReaK

FReaK has been published as a tool paper at the 27th ACM International Conference on Hybrid Systems: Computation and Control (HSCC) 2024, available [here](https://dl.acm.org/doi/abs/10.1145/3641513.3650141).

If you cite FReak, please cite

Bak, S., Bogomolov, S., Hekal, A., Kochdumper, N., Lew, E., Mata, A. and Rahmati, A., 2024, May. Falsification using Reachability of Surrogate Koopman Models. In Proceedings of the 27th ACM International Conference on Hybrid Systems: Computation and Control (pp. 1-13)

Bibtex:
```
@inproceedings{bak2024falsification,
  title={Falsification using Reachability of Surrogate Koopman Models},
  author={Bak, Stanley and Bogomolov, Sergiy and Hekal, Abdelrahman and Kochdumper, Niklas and Lew, Ethan and Mata, Andrew and Rahmati, Amir},
  booktitle={Proceedings of the 27th ACM International Conference on Hybrid Systems: Computation and Control},
  pages={1--13},
  year={2024}
}
```

## Contact

