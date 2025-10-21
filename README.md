# ScalingFromSkypackage

<!-- badges: start -->

[![R-CMD-check](https://github.com/ForestScaling/ScalingFromSkypackage/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ForestScaling/ScalingFromSkypackage/actions/workflows/R-CMD-check.yaml)

<!-- badges: end -->


`ScalingFromSkypackage` is an R package designed to estimate forest size-abundance distributions using remote sensing imagery and field survey data. It provides a workflow for integrating structural forest measurements, canopy height models, and statistical modeling to produce estimates of tree abundance and size distribution parameters.  

The package implements methods for estimating the Pareto α parameter, which characterizes the negative power-law relationship commonly observed in forest tree size distributions, and for predicting total tree abundance while accounting for factors such as canopy structure and sampling corrections.  

## Key Features

- Flexible Bayesian modeling of size-abundance distributions  
- Estimation of total tree counts with LAI and breakpoint adjustments  
- Compatibility with remote sensing and field-collected forest plot data  
- Example dataset for testing and demonstration purposes  

## Included Data

- **`harv_data_sample`**: A curated sample dataset derived from the Harvard Forest CTFS-ForestGEO plot (HF253) and NEON Airborne Observation Platform (AOP). This dataset includes tree crown coordinates, canopy structural attributes, derived DBH estimates, and metadata. It is intended for demonstration, testing, and educational purposes.  

## Purpose

This package makes previously published methods for estimating forest size-abundance distributions more accessible to researchers and practitioners. It is particularly useful for:

- Forest ecologists analyzing structural changes across plots or regions  
- Researchers integrating remote sensing data with field measurements  
- Teaching and demonstration of forest structural modeling workflows  

## License

MIT License
