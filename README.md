# Thesis Chapter 3
Thesis Chapter Code:

This repository contains the code used in chapter 3 of the thesis: Understanding how mutability facilitates survival of alternating selection and botlenecks by the major food-borne pathogen Campylobacter jejuni (2023)
and from:
Jonathan Holmes, Lickson Munjoma, Christopher D. Bayliss,
Novel method for prediction of combinatorial phase-variable gene expression states,
MethodsX,
Volume 11,
2023

The code used is listed as follows:

Algorithm_V1.1.py            <- Sweep corection algorithm

InSilicoAnalysis.py          <- In silico comparison of single colonies and sweep data collection\n

In_Vitro_Data.csv            <- In vitro experimental data  

in_vitro_data_comp.py        <- In vitro comparison script for sweep and single colonies

missing_data_in_invitro.py   <- Simulation of sweep errors for the in vitro data


The sweep - correction algorithm can be run using:

python3 Algorithm_V1.1.py input_file.csv

where input_file.csv is of identical format to that of In_Vitro_Data.csv  

the remaining scripts can be called as standalone scripts (requiring the Algorithm_V1.1.py and In_Vtiro_Data.py within the same directory) with parameters edited within the script text file.

The libraries required for the running of each script can be found at the top of each .py document and comments annotated where changes can be made to the runnning of the script.  



