# Matlab Coral Bleaching Model

This contains a Matlab version of the coral bleaching model contained within the C-based CSIRO Environmental Modelling Suite [(GitHub - csiro-coasts/EMS: Environmental Modelling Suite)](https://github.com/csiro-coasts/EMS). EMS has been used for modelling circulation, sediments, chemistry and ecology on coral reef systems (eReefs Project, ereefs.info). It contains a mechanistic model of the coral-symbiont relationship that considers temperature-mediated build-up of reactive oxygen species due to excess light, leading to zooxanthellae expulsion. The model explicitly represents the coral host biomass, as well as zooxanthellae biomass, intracellular pigment concentration, nutrient status, and the state of reaction centres and the xanthophyll cycle. Photophysiological processes represented include photoadaptation, xanthophyll cycle dynamics, and reaction centre state transitions.

The original Matlab version was written by Mark Baird and Julie Terp Jørgensen.

A version of the model was developed for a coral tank laboratory experiment. See Bleach_SCU_lab_exp_2024.m and the configuration to run this (Run_Bleach_SCU_lab_exp_2024.m) for details on alterations made. 

References:

Baird, M. E., K. Wild-Allen, J. Parslow, M. Mongin, B. Robson, J. Skerratt, F. Rizwi, M. Soja-Woźniak, E. Jones, M. Herzfeld, N. Margvelashvili, J. Andrewartha, C. Langlais, M. Adams, N. Cherukuru, S. Hadley, P. Ralph, T. Schroeder, A. Steven, U. Rosebrock, L. Laiolo, M. Gustafsson, and D. Harrison (2020). CSIRO Environmental Modelling Suite (EMS): Scientific description of the optical and biogeochemical models (vB3p0). Geoscientific Model Development.13:4503-4553.

Baird, M. E., M. Mongin, F. Rizwi, L. K. Bay, N. E. Cantin, M. Soja-Wozniak and J. Skerratt (2018) A mechanistic model of coral bleaching due to temperature-mediated light-driven reactive oxygen build-up in zooxanthellae. Ecol. Model 386: 20-37.
