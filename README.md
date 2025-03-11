# Improved effectiveness of vaccination campaigns against rabies by reducing spatial heterogeneity in coverage

Authors: Elaine A Ferguson, Ahmed Lugelo, Anna Czupryna, Danni Anderson, Felix Lankester, Lwitiko Sikana, Jonathan Dushoff, Katie Hampson

This repository includes the code and de-identified data used in Ferguson et al. 2025 to analyse of the impacts of rabies vaccination in Serengeti District over 2002-2022. 

The script workflow.R sources the analysis scripts in the correct order, with annotations giving a brief description of the purpose of each.

To de-identify the spatial data on rabies transmission (which include locations of households of dog owners or persons bitten by rabid animals), x and y coordinates were jittered uniformly within 1km. As a result, distance kernels, inferred incursions, and village-level models estimated using the analysis script in this repository will differ slightly from those presented in the paper. We have therefore included deidentified results produced using the exact locations in the output folder, which are read in at the appropriate points to produce the plots etc. as published.

For more information: elaine.ferguson@glasgow.ac.uk
