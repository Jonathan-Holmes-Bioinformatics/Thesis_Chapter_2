# -----------------------------------------------------------------------------------------------
# 
#
#
# Jonathan Holmes
#
# Permission to use, copy, modify, and/or distribute this software or any part thereof for any
# purpose with or without fee is hereby granted provided that:
#     (1) the original author is credited appropriately in the source code
#         and any accompanying documentation
# and (2) that this requirement is included in any redistribution.
#
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS
# SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL
# THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT,
# NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE
# OF THIS SOFTWARE.
#
# E-mail: jh795@leicester.ac.uk
# ------------------------------------------------------------------------------------------------
# 
#
# This script should rely on no libraries outside of the default python libraries and as such should
# require no additional installation. To run the programme navigate to the direcotry containing your
# input file then run the programme as follows:
#
# 			python3 in_vitro_data_comp.py
#
# 
#
# ------------------------------------------------------------------------------------------------


### Required functions - installed through conda in orginal code

import os
import itertools
import parsimony_v2
import random
import sklearn
from scipy.stats import entropy
from sklearn import metrics
import scipy
from scipy.spatial import distance



## in built functions

def phasome(A):
	E = []
	for B in A:
		C = ""
		for D in B:
			C = C + str(D)
		E.append(C)
	return E

def to_numerical(A):
	B = []
	for i in A:	
		B.append(i/100)
	return B

def normalise_to_sum(A):
	maximum = sum(A)
	B = []
	for i in A:
		B.append(i/maximum)
	return B

def to_float(A):
	B = []
	for i in A:
		B.append(float(i))
	return B


def to_string(A):
	B = []
	for i in A:
		B.append(str(i))
	return B

def sum_of_squares(A,B):

	value = 0
	
	i = 0
	for I in A:
		value = value + (float(I) - float(B[i]))*(float(I) - float(B[i]))
		i = i + 1

	return value




## Change input file here
Input_file = open("input_file").read().split("\n")

## Defines the Phasotype distribution and calculates the sweep for that data.
Phasotypes = phasome(list(itertools.product([0, 1], repeat=len(Input_file[0].split(",")[1:]))))
Distribution = [0]*len(Phasotypes)
for line in Input_file[2:]:
	l = line.split(",")[1:]
	pt = "".join(l)
	if pt in Phasotypes:
		Distribution[Phasotypes.index(pt)] = Distribution[Phasotypes.index(pt)] + 1
Distribution_in_vitro = normalise_to_sum(Distribution)
Sweep = to_numerical(to_float(Input_file[1].split(",")[1:]))


## Output_file
Output = open("Corrected_Divergences.csv","w")
## Contains output_data
Matrix = []


## The range of single colony samples to be compared so: for sample in range(1,10) will test single colonies sample sizes of 1,2,3,4,5,6,7,8,9 and 10
for sample in range(1,11):
	print(sample)
	Line = [sample]
	
	## iter refers to the number of replciates for each sample size so: for iter in range(0,20): will run 20 iterations (0,30) will run 30 iterations etc. 
	samples = []
	for iter in range(0,20):
		samples.append(random.choices(Phasotypes, weights=Distribution_in_vitro, k=sample))
	## Samples are corrected using sweep data		
	corrected_pop = []
	sampled_pop = []
	for d in samples:
		distribution = [0]*len(Phasotypes)
		for D in d:
			distribution[Phasotypes.index(D)] = distribution[Phasotypes.index(D)] + 1
		sampled_pop.append(distribution)
		corrected_pop.append(parsimony_v2.parsimony_algorithm(normalise_to_sum(distribution), Sweep, Phasotypes)[1][-1])

	
	
	sampled_div = []
	correct_div = []

	#### Divergence from "true" values estimated 
	for i in range(0,len(sampled_pop)):

		sampled_div.append(distance.jensenshannon(normalise_to_sum(sampled_pop[i]), Distribution_in_vitro))
		correct_div.append(distance.jensenshannon(corrected_pop[i], Distribution_in_vitro))
	Line = Line + sampled_div + correct_div
	Matrix.append(Line)
	

## Data written to output file. 
for l in Matrix:
	Output.write(",".join(to_string(l)) + "\n")


Output.close()


































 



































