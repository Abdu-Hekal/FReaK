# Koopman Falsification

## Overview

A MATLAB toolkit for falsification using Koopman surrogate models. 

## Table of Contents

- [Introduction](#introduction)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Usage](#usage)
- [Publications](#publications)
- [Contact](#contact)

## Introduction

[Provide a brief introduction to the Koopman Falsification project. Include its purpose, goals, and any relevant background information.]

## Prerequisites

Koopman falsification requires the following prerequisites:

- MATLAB - You can download MATLAB from [MathWorks website](https://www.mathworks.com/).
- Breach - Download Breach toolbox from github repository [here](https://github.com/decyphir/breach/tree/master)
- CORA(v2022) - Download CORA v2022 from the archive page on the official [CORA website](https://tumcps.github.io/CORA/). Follow the instructions in the corresponding manual for installation.
- YALMIP - See CORA2022 manual for best installation description
- Python(>=3.8) - Install Python version 3.8 or higher.
- AutoKoopman - Install AutoKoopman directly from PyPI using: `pip install autokoopman`
- Gurobi Optimization Toolbox - We use Gurobi solver for backend optimization, though other solvers may be used. Follow the instructions provided on the [Gurobi website](https://support.gurobi.com/hc/en-us/articles/4533938303505-How-do-I-install-Gurobi-for-Matlab-) to download and install Gurobi for matlab

## Installation

We provide a setup script to help with installation. The script automatically installs prerequisites where possible 
or outlines which are missing via the command window. The script also makes necassary modifications to external toolboxes
and removes any conflicting files from path. Add the Koopman falsification folder and subfolders to path and run setup file
**setupKF.m**

```matlab
% Run setup file
setupKF
```

Note that the script also ensure that only necessary prerequisites for external toolboxes (like CORA) are installed.
For full functionality of each external toolbox, refer to their corresponding instructions.

## Usage

## Publications

## Contact

