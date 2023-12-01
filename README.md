# BNPower - a power calculation tool for data-driven network analysis for whole-brain connectome data

## Overview
BNPower is a specialized tool designed for power calculations in data-driven network analysis of whole-brain connectome data. This tool accompanies our paper, "BNPower - A Power Calculation Tool for Data-Driven Network Analysis in Whole-Brain Connectome Data," published in [XXX journal].

## Features of BNPower
BNPower is equipped with two distinct tabs catering to different research designs:
1. T-test Tab: Ideal for studies requiring a two-sample test, such as case versus control comparisons. This tab facilitates the comparison of two distinct groups within your connectome data.
2. Regression Tab: Suited for studies where the primary variable of interest is continuous, like age. This tab supports the examination of linear relationships between continuous variables and network features in connectome data.

<img src="https://github.com/bichuan0419/brain_connectome_power_tool/assets/43563121/d87fd1e2-8898-4fa5-b686-33f753dfbc77" width="400">
<img src="https://github.com/bichuan0419/brain_connectome_power_tool/assets/43563121/17034740-c9cd-4dbd-855d-79bbe0dee084" width="400">


## How to use the BNPower
BNPower is developed in MATLAB and offers multiple usage options depending on whether MATLAB is installed on your system. Follow these steps to get started:
1. For Users Without MATLAB:
Download and run the Webinstaller file from BNPower/for_redistribution/MyAppInstaller_web.exe. This will install the MATLAB Runtime required to run BNPower.
After installing the MATLAB Runtime, use BNPower/for_redistribution_files_only/BNPower.exe to launch the application.
2. For Users With MATLAB Installed:
You can directly run BNPower/for_redistribution_files_only/BNPower.exe.
3. Alternatively, download and unzip the complete BNPower package and use BNPower.mlapp within MATLAB.

User can simply click "Run" button for power calculation. A "Show Example Network" button is available for users to inspect the sample inference matrix. A concise video tutorial for BNPower is available at https://youtu.be/Hd1splkxrhU

## Inputs
There are 3 main categories of the input variables: 1. the input that governs the graph structure; 2. The input that are requried for classical univariate power calculation (e.g., sample size, effect size, alpha level); 3. The inputs that affects power in the simulation-based power calculation process (e.g., number of datasets used, number of permutation test for each dataset).

### The pipeline of BNPower is shown here:
<img src="https://github.com/bichuan0419/brain_connectome_power_tool/assets/43563121/5a2f86c1-f2e4-4628-bb2f-b42b26268508" width="600">

### A summary of the inputs can be found below:
<img src="https://github.com/bichuan0419/brain_connectome_power_tool/assets/43563121/f75f972b-a703-4636-97a2-6383ecc59fcc" width="600">

### The network-level statistical analysis procedure can be found here
<img src="https://github.com/bichuan0419/brain_connectome_power_tool/assets/43563121/abfcc4a8-045e-4c3d-8c15-808b612b9817" width="600">

## Outputs
Two outputs will be recorded for BNPower: The network-level statistical power with 95% confidence interval (as indicated in the bolded, red textfield), and the time elapsed to calculate the power.


## Compatibility
The software is compatible with MATLAB 2019b and newer versions. The required toolboxes are shown below:
<img src="https://github.com/bichuan0419/brain_connectome_power_tool/assets/43563121/f490674e-d8ff-4c38-b87c-9d95a86f1d96" width="600">

## Runtime
The estimated runtime for BNpower varies depending on the specific input settings. We have compiled a table for the expected runtime of BNPower using a Windows PC equipped with a 13th Gen Intel(R) Core(TM) i5-1340P processor (1.90 GHz) and 32.0 GB RAM (31.6 GB usable). The system operated on a 64-bit Windows OS with an x64-based processor. 

<img src="https://github.com/bichuan0419/brain_connectome_power_tool/assets/43563121/ac47b554-a349-40f5-b695-6508c13dbbdc" width="600">


## License
BNpower is developed under the GNU Licesnce.
