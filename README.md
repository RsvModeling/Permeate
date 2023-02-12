# Permeate
PERMutation basEd Analysis of mulTiple Endpoints

This package implements the methods described in the manuscript: "Incorporating data from multiple endpoints in the analysis of clinical trials: example from RSV vaccines",
available at: https://www.medrxiv.org/content/10.1101/2023.02.07.23285596v2.

## Abstract: 
To achieve licensure, interventions typically must demonstrate efficacy against a primary outcome in a randomized clinical trial. However, selecting a single primary outcome a priori is challenging. Incorporating data from multiple and related outcomes might help to increase statistical power in clinical trials. Inspired by real-world clinical trials of interventions against respiratory syncytial virus (RSV), we examined methods for analyzing data on multiple endpoints.

We simulated data from three different populations in which the efficacy of the intervention and the correlation among outcomes varied. We developed a novel permutation-based approach that represents a weighted average of individual outcome test statistics (varP) to evaluate intervention efficacy in a multiple endpoint analysis. We compared the power and type I error rate of this approach to two alternative methods: the Bonferroni correction (bonfT) and another permutation-based approach that uses the minimum P-value across all test statistics (minP).

When the vaccine efficacy against different outcomes was similar, VarP yielded higher power than bonfT and minP; in some scenarios the improvement in power was substantial. In settings where vaccine efficacy was notably larger against one endpoint compared to the others, all three methods had similar power.

Analyzing multiple endpoints using a weighted permutation method can increase power while controlling the type I error rate in settings where outcomes share similar characteristics, like RSV outcomes. We developed an R package, PERMEATE, to guide selection of the most appropriate method for analyzing multiple endpoints in clinical trials.

## How to install the PERMEATE package
To install the package from github:
install.packages("devtools")
devtools::install_github("weinbergerlab/Permeate")

# Instructions
The user can either run one of the scenarios described in the manuscript, or can input their own parameters to run the simulation framework. 
## Main parameters:
- N.sim: number of simulated datasets
- N.outcomes: number of outcomes (minimum number is 2)
- RR: risk ratio for each outcome
- prop.outcome: incidence of each outcome
- N1: number of individuals in the control arm
- N2: number of individuals in the treatment arm
## Pre-specified scenarios: 
1. "setting1": scenario A) in the manuscript. Parameters are set as: N.sim=1000,N.outcomes=3,RR=c(0.60,0.6,0.7),prop.outcome=c(.22,.20,.12),N1=200,N2=200.
2. "setting2": scenario B) in the manuscript. Parameters are set as: N.sim=1000,N.outcomes=3,RR=c(0.25,0.4,0.6),prop.outcome=c(.05,.02,.03),N1=496,N2=994.
2. "setting3": scenario C) in the manuscript. Parameters are set as: N.sim=1000,N.outcomes=3,RR=c(0.60,0.55,0.5),prop.outcome=c(.02,.04,.01),N1=1430,N2=2765.
## How to run the package: 
Step 1: The first function to run is the set_setting function, example:set_setting("setting1"). If "custom", the users needs to input all the parameters described in the previous section. Example: set_setting("custom",N.sim=1000,N.outcomes=3,RR=c(0.60,0.55,0.5),prop.outcome=c(.02,.04,.01),N1=1430,N2=2765)

Step 2: Run the main function. The inputs are: the setting specified in Step 1, the setting name, correlation among the outcomes, and directory where to store the results. Example: main_run(setting,setting_name,corr=cor_l[i],dir)

Please refer to the vignette for detailed example.


