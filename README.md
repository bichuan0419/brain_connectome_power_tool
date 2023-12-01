# BNPower - a power calculation tool for data-driven network analysis for whole-brain connectome data
This code accompanies the paper "BNPower - a power calculation tool for data-driven network analysis for whole-brain connectome data" on XXX journal
## Description of BNPower
There are two tabs within BNPower, researchers can use "T-test" in the study designs where a two-sample test (case vs control) is needed, users can also use "Regression" tab for the studies where the predictor-of-interest is a continuous variable (e.g. age).

[PASTE FIGURE OF THE GUI HERE]

## INPUTS
There are 3 main categories of the input variables: 1. the input that governs the graph structure; 2. The input that are requried for classical univariate power calculation (e.g., sample size, effect size, alpha level); 3. The inputs that affects power in the simulation-based power calculation process (e.g., number of datasets used, number of permutation test for each dataset).

### The pipeline of BNPower is shown here:
[PASTE figure]

### A summary of the inputs can be found below:
[PASTE FIGURE]

### The network-level statistical analysis procedure can be found here
[PASTE FIGURE]

## How to use the tool
BNPower is developed and compiled in MATLAB, so to use the tool, users can either
1. Use the Webinstaller file 
