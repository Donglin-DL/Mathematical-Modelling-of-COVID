# Mathematical-Modelling-of-COVID
This is the code for Paper "Accurate Estimation of the Infected Proportion during the
First Wave of COVID-19 in the Province of Alberta, Canada:
A Mathematical Modeling Study."

Data used for this paper:

Daily covid case data, daily covid testing data, and daily covid death data. We thank the Analytics and Performance Reporting Branch at Alberta Health for data support.

To run the code and view the results, please follow these steps:
1. Run death_data.m
2. Run ModelScript.m : You will be fitting the mathematical model and obtain model outputs such as posterior distribution of each parameters,the estimated susceptibles, the estimated infected population, and how well the model fits the cases.
3. Run Results_out.m : You will be getting more information about the model such as the proportions of number of hidden infections and newly diagnosed infection, sensitivity analysis of certain outcomes.

