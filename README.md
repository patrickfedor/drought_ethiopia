# Ethiopia Climate Modeling Analysis

## Project Description:

This project analyzes the probability of extreme drought across Ethiopia
using climate models. The analysis relies on an ensemble of climate
models with a baseline reference period from 1985-2015. The focus period
for this project is from 2041 to 2080. This project is meant to provide
insights into Ethiopia’s future climate vulnerabilities.

# Methodology and Data

The ensemble includes 18 climate models, with two models excluded due to
incomplete data. To estimate evapotranspiration, Hargreaves formulation
for Potential Evapotranspiration (PET) was applied. This method uses
temperature and precipitation data to offer a reliable estimate in areas
where data availability is limited. Water balance calculations were as a
rolling sum and performed with a 12-month integration window.

Extreme drought in the context of this project is defined as an event
with a historical probability of 10% during the baseline period. This is
consistent with methodologies like the Standardized Precipitation
Evapotranspiration Index (SPEI) by Serrano et al. Unlike the traditional
approach that fit statistical distributions, this project calculated
deviations directly using quantiles, providing an easily interpretable
measure of changing drought probabilities.

![Probability of Extreme Drought](figures/prob_extr_drought_plot.png)
