R Code for "Profile monitoring via a nonparametric Bayesian online learning framework"
-------------------------------------------------------------------------------------------------------------
by Daewon Yang
-------------------------------------------------------------------------------------------------------------
Function.R 
 -> functions to implement the proposed model

DataGen.R 
-> Generate the simulation data (Scenario 1, 3, 4)

DataGen__sc2.R 
-> Generate the simulation data (Scenario 2)

Part0___Choose_Control_limit.R
-> Code to conduct Monte Carlo simulatons for Phase I - proposed model (Scenario 1, 3, 4)

Part0___Choose_Control_limit__sc2.R
-> Code to conduct Monte Carlo simulatons for Phase I - proposed model (Scenario 2)

Part0___Choose_Control_limit__Freq.R
-> Code to conduct Monte Carlo simulatons for Phase I - competing model using K-means algorithm

Part1___Result1_for_Proposed_model.R
-> Code for Phase II monitoring - proposed model 

Part2__Result2_3_Clustering_ARL1.R
-> Code to conduct Monte Carlo simulation for Phase II - proposed model 

Part2__Result2_3_Clustering_ARL1__Freq.R
-> Code to conduct Monte Carlo simulation for Phase II - competing model using K-means algorithm

CalculateL.R
-> Code to calculate an approporiate L with ARL0 = 100 

CalculateARL1.R
-> Code to calculate ARL1 using the Monte Carlo simulation result 
-------------------------------------------------------------------------------------------------------------