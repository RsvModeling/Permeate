# Permeate
PERMutation basEd Analysis of mulTiple Endpoints

This package implements the methods described in the manuscript: "Incorporating data from multiple endpoints in the analysis of clinical trials: example from RSV vaccines",
available at:.
The user can either run one of the scenarios described in the manuscript, or can input their own parameters to run the simulation framework. 
## Main parameters:
- N.sim: number of simulated datasets
- N.outcomes: number of outcomes (minimum number is 2)
- RR: risk ratio for each outcome
- prop.outcome: incidence of each outcome
- N1: number of individuals in the control arm
- N2: number of individuals in the treatment arm
## Pre-specified scenarios: 
1. "setting1": this scenario corresponds to scenario A) in the manuscript. Parameters are set as: N.sim=1000,N.outcomes=3,RR=c(0.60,0.6,0.7),prop.outcome=c(.22,.20,.12),N1=200,N2=200.
2. "setting2": this scenario corresponds to scenario B) in the manuscript. Parameters are set as: N.sim=1000,N.outcomes=3,RR=c(0.25,0.4,0.6),prop.outcome=c(.05,.02,.03),N1=496,N2=994.
2. "setting3": this scenario corresponds to scenario C) in the manuscript. Parameters are set as: N.sim=1000,N.outcomes=3,RR=c(0.60,0.55,0.5),prop.outcome=c(.02,.04,.01),N1=1430,N2=2765.
## How to run the package: 

Step 1: The first function to run is the set_setting function. If the user wants to run one of the scenarios described in the manuscript, for example setting1, it is sufficient to run: set_setting("setting1"). 
Otherwise, the user must specify "custom" and all the parameters described in the previous section. Example: set_setting("custom",N.sim=1000,N.outcomes=3,RR=c(0.60,0.55,0.5),prop.outcome=c(.02,.04,.01),N1=1430,N2=2765)-
