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
You can directly run BNPower/for_redistribution_files_only/BNPower.exe. Alternatively, download and unzip the complete BNPower package and use BNPower.mlapp within MATLAB.

### Operating BNPower
* To perform a power calculation, simply click the "Run" button in the application.
* Use the "Show Example Network" button to view a sample inference matrix, which can help familiarize you with the tool’s functionality.
* A concise video tutorial for BNPower is available at this YouTube link for additional guidance.

## Inputs for BNPower

BNPower requires three main categories of input variables:

1. **Graph Structure Inputs**: These inputs determine the structure of the graph in the analysis. They define how the network nodes (N) and edges are organized and interact within the connectome data.

2. **Classical Univariate Power Calculation Inputs**: This category includes the traditional parameters necessary for univariate power calculations, such as:
   - Sample Size: The number of observations or data points in each group.
   - Effect Size: The anticipated size of the effect or difference you are trying to detect.
   - Alpha Level: The significance threshold, typically set at 0.05, which determines the probability of a Type I error (false positive).

3. **Simulation-Based Power Calculation Inputs**: These inputs are crucial for the simulation aspect of power calculations and include:
   - Number of Datasets Used: The quantity of datasets utilized in the simulation process.
   - Number of Permutation Tests per Dataset: The frequency of permutation tests conducted for each dataset to assess the power accurately.

### A summary of the inputs:
<img src="https://github.com/bichuan0419/brain_connectome_power_tool/assets/43563121/f75f972b-a703-4636-97a2-6383ecc59fcc" width="600">

### The pipeline of BNPower:
<img src="https://github.com/bichuan0419/brain_connectome_power_tool/assets/43563121/5a2f86c1-f2e4-4628-bb2f-b42b26268508" width="600">

### The network-level statistical analysis procedure:
<img src="https://github.com/bichuan0419/brain_connectome_power_tool/assets/43563121/abfcc4a8-045e-4c3d-8c15-808b612b9817" width="600">

## Outputs of BNPower

BNPower provides two primary outputs:

1. **Network-Level Statistical Power with 95% Confidence Interval**: This is displayed in a bolded, red text field within the application. It indicates the statistical power of your network analysis along with its 95% confidence interval, offering a clear understanding of the robustness of your results.

2. **Time Elapsed for Power Calculation**: The application also records and displays the time taken to calculate the power, providing insight into the efficiency of the analysis process.



## Compatibility
The software is compatible with MATLAB 2019b and newer versions. The required toolboxes are shown below:
<img src="https://github.com/bichuan0419/brain_connectome_power_tool/assets/43563121/f490674e-d8ff-4c38-b87c-9d95a86f1d96" width="600">

## Runtime
The estimated runtime for BNpower varies depending on the specific input settings. We have compiled a table for the expected runtime of BNPower using a Windows PC equipped with a 13th Gen Intel(R) Core(TM) i5-1340P processor (1.90 GHz) and 32.0 GB RAM (31.6 GB usable). The system operated on a 64-bit Windows OS with an x64-based processor. 

<img src="https://github.com/bichuan0419/brain_connectome_power_tool/assets/43563121/ac47b554-a349-40f5-b695-6508c13dbbdc" width="600">


## License

BNPower is developed under the GNU General Public License. This open-source license ensures that users have the freedom to run, study, share, and modify the software. For more details, please refer to the [GNU General Public License](https://www.gnu.org/licenses/gpl-3.0.en.html).

