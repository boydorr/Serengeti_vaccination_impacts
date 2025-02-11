# Reducing spatial heterogeneity in coverage improves the effectiveness of dog vaccination campaigns against rabies

https://www.biorxiv.org/content/10.1101/2024.10.03.616420v3

Authors: Elaine A Ferguson, Ahmed Lugelo, Anna Czupryna, Danni Anderson, Felix Lankester, Lwitiko Sikana, Jonathan Dushoff, Katie Hampson

This repository includes the code and de-identified data used in this preprint to analyse of the impacts of rabies vaccination in Serengeti District, 2002-2022. The script workflow.R sources the analysis scripts in the correct order, with a brief description of the purpose of each.

To de-identify the spatial data on rabies transmission (which include locations of households of dog owners or persons bitten by rabid animals), x and y coordinates were jittered uniformly within 1km. As a result, distance kernels, inferred incursions, and village-level models from this repository will differ slightly from those presented in the paper.

For more information: elaine.ferguson@glasgow.ac.uk
