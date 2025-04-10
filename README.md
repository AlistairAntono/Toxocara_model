# Toxocara_model
This Repository Contains the Code for Compartmental Model of Zoonotic Transmission of Toxocara

The model code is designed to be run with specifically formatted input files. These can be found in the repository using input data which has been used in the paper "Development of a Dynamic Stochastic Compartmental Model of Zoonotic Toxocariasis Transmission" submitted to Zoonoses and Public Health, 2025.

The model requires Input_demographics, Input_parameters, and transmission_params dataframes to run correctly. 

Input_demographics (see associated file in the repository for formatting) contains demographic information for the country(ies) or population(s) to be used for the analysis divded by age group. Age groups are in 10 year intervals (0-9; 10-19; 20-29; etc) up to 80+.

Input_parameters (see associated file in the repository for formatting) contains information for the country(ies) or population(s) to be analysed: the estimated dog and cat populations for the preceding 3 years (used to estimate the population growth rate within the model code), the initial human population, the yearly human population growth rate, the prevalence of toxocara in dogs and cats, the initial proportion of the human population to be placed in the recovered compartment (expressed as a decimal - corresponds to the seroprevalence estimate for that country/population), the estimated number of infections (lower, mean, upper CI), and economic cost estimates for each manifestation of toxocariasis (visceral, ocular, common), and the region the country is in. This is necessary to assign the regional transmission rate.

transmission_params contains regional transmission rates for toxocariasis derived from iterative model fitting using the limited memory Broyden-Fletcher-Goldfarb-Shanno method with bounds (L-BFGS-B) using the optimx function in R. Code used to estimate the transmission rates is also provided separately in this repository, and can be run using the same input files (minus the transmission_params files)

The model runs for 1000 iterations, and collects data every month for 120 months of simulation.

The model then outputs population data, summary statistics, and a series of figures which allow visualisation of the results, which are saved to the working directory.

Within the monthly simulation loop, the model samples the transmission rate values across a PERT distribution for each month of the simulation. The model then stores the sampled values for the iteration, which can be viewed in the output files.

Exposures are then calculated based on the monthly transmission rate, and a scaling factor is applied based on the age breakdown of the population - to reflect a higher likelihood of transmission within certain age groups based on UK clinical data. Exposures by age group are then summed.

The symptomatic proportion of exposures is then calculated by sampling across a PERT distribution for the symptomatic cases range (20-40%, mid 30%). Total new exposures is then multiplied by the sampled % for the month. Severe cases are calculated based on the severe cases percentage (fixed parameter)

The populations compartments are then updated accordingly, and the calculated monthly growth rate is applied.

Monthly economic costs are obtained by the multiplication of the estimated monthly costs by the number of cases.

This is done for each country within the Input files.
