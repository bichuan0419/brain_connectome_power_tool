# BNPower - a power calculation tool for data-driven network analysis for whole-brain connectome data
This code accompanies the paper "BNPower - a power calculation tool for data-driven network analysis for whole-brain connectome data" on XXX journal
## Description of BNPower
There are two tabs within BNPower, researchers can use "T-test" in the study designs where a two-sample test (case vs control) is needed, users can also use "Regression" tab for the studies where the predictor-of-interest is a continuous variable (e.g. age).

![Fig4a](https://github.com/bichuan0419/brain_connectome_power_tool/assets/43563121/d87fd1e2-8898-4fa5-b686-33f753dfbc77)

![Fig4b](https://github.com/bichuan0419/brain_connectome_power_tool/assets/43563121/17034740-c9cd-4dbd-855d-79bbe0dee084)

## How to use the tool
BNPower is developed and compiled in MATLAB, so to use the tool, users can either
1. Use the Webinstaller file in BNPower/for_redistribution/MyAppInstaller_web.exe to install Matlab Runtime, if the user does not have Matlab installed.
2. If the user has Matlab installed, he/she can use /BNPower/for_redistribution_files_only/BNPower.exe or
3. If the user has Matlab installed, he/she can use BNPower.mlapp after downloading/unzipping the complete package on Matlab.

User can simply click "Run" button for power calculation. A "Show Example Network" button is available for users to inspect the sample inference matrix. A concise video tutorial for BNPower is available at https://youtu.be/Hd1splkxrhU

## INPUTS
There are 3 main categories of the input variables: 1. the input that governs the graph structure; 2. The input that are requried for classical univariate power calculation (e.g., sample size, effect size, alpha level); 3. The inputs that affects power in the simulation-based power calculation process (e.g., number of datasets used, number of permutation test for each dataset).

### The pipeline of BNPower is shown here:
![Fig3a](https://github.com/bichuan0419/brain_connectome_power_tool/assets/43563121/5a2f86c1-f2e4-4628-bb2f-b42b26268508)

### A summary of the inputs can be found below:
![Fig3b](https://github.com/bichuan0419/brain_connectome_power_tool/assets/43563121/f75f972b-a703-4636-97a2-6383ecc59fcc)


### The network-level statistical analysis procedure can be found here
![Fig2](https://github.com/bichuan0419/brain_connectome_power_tool/assets/43563121/abfcc4a8-045e-4c3d-8c15-808b612b9817)


## Compatibility
The software is compatible with MATLAB 2019b and newer versions. The required toolboxes are shown below:
![Required toolboxes](https://github.com/bichuan0419/brain_connectome_power_tool/assets/43563121/f490674e-d8ff-4c38-b87c-9d95a86f1d96)

## Runtime
The estimated runtime for BNpower varies depending on the specific input settings. We have compiled a table for the expected runtime of BNPower using a Windows PC equipped with a 13th Gen Intel(R) Core(TM) i5-1340P processor (1.90 GHz) and 32.0 GB RAM (31.6 GB usable). The system operated on a 64-bit Windows OS with an x64-based processor. 
![image](https://github.com/bichuan0419/brain_connectome_power_tool/assets/43563121/ac47b554-a349-40f5-b695-6508c13dbbdc)


## License
BNpower is developed under the GNU Licesnce.
