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
# 			python3 Sweep_Error.py
#
# 
#
# ------------------------------------------------------------------------------------------------


### Required functions - installed through conda in orginal code

import os
import sys
from scipy.spatial import distance
import random
import itertools


# Input = input file
# index = formatting
# Size = number of single colonies for sweep correction
# Error Range  = range of sweep errors to be tested
# repeats = number of replicates

Input = open("Example_Input.csv").read().split("\n")
index = Input[0:1]
Size = 7
Error_range = [0,1,2,3,4,5,10,20,50,100]  #out of 100
repeats = 20
True_Distribution = []


## Inbuilt functions

def to_float(A):
	B = []
	for i in A:
		B.append(float(i))
	return B

def normalise_to_sum(A):
	maximum = sum(A)
	B = []
	for i in A:
		B.append(i/maximum)
	return B

def phasome(A):
	E = []
	for B in A:
		C = ""
		for D in B:
			C = C + str(D)
		E.append(C)
	return E

def to_string(A):	
	B = []
	for i in A:
		B.append(str(i))
	return B


# Building phasotype distribution 
Phasotypes = phasome(list(itertools.product([0, 1], repeat=len(Input[0].split(",")[1:]))))

Distribution = [0]*len(Phasotypes)
for line in Input[2:]:
	l = line.split(",")[1:]
	pt = "".join(l)
	if pt in Phasotypes:
		Distribution[Phasotypes.index(pt)] = Distribution[Phasotypes.index(pt)] + 1

Distribution_in_vitro_1 = normalise_to_sum(Distribution)
Distribution_in_vitro = []
for i in Distribution_in_vitro_1:
	if i == 0.0:
		Distribution_in_vitro .append(0.000000001)
	else:
		Distribution_in_vitro .append(i)



# Builds a baseline of single colonies to compare against - set to 16 - to change this alter for iter in range(0,16): to for iter in range(0,N): where N is the number of single colonies to be tested
baseline = []
for i in range(0,20):
	I = []
	for iter in range(0,16):
		I.append(random.choice(Input[3:-2]))
	baseline.append(I)
phasotype_baseline = []

for s in baseline:
	base_pv = [0]*len(Phasotypes)
	for S in s:
		phase = "".join(S.split(",")[1:])
		base_pv[Phasotypes.index(phase)] = base_pv[Phasotypes.index(phase)] + 1
	
	phasotype_baseline.append(normalise_to_sum(base_pv))

colony_only = []

for col_pv in phasotype_baseline:
	colony_only.append(str(distance.jensenshannon(col_pv,Distribution_in_vitro)))



## Outputs the error adjusted data from the error range.


outputs =  []


for error in Error_range:
	matrix = []
	for rep in range(0,repeats):
		samples = []
		for iter in range(0,Size):
			samples.append(random.choice(Input[3:-2]))

		matrix.append(samples)
	## makes random choice on adding or removing from sweep values
	sweeps = to_float(Input[1].split(",")[1:])

	choice = [error, 0 - error]

	# Error adjusted sweeeps 
	sweep_matrix = []
	for rep in range(0,repeats):
		Error_sweeps =[]
		for sweep in sweeps:
			alter = sweep + random.choice(choice)
			if 0 < alter and alter < 100:
				Error_sweeps.append(alter)
			if 100 <= alter:
				Error_sweeps.append(100)
			if alter <= 0:
				Error_sweeps.append(0)		
		sweep_matrix.append(Error_sweeps)		

	output_error = []

	# Runs correction on adjusted sweeps
	for run in range(0,len(matrix)):
		temp = open("temp_input.csv","w")
		temp.write(index[0] + "\n")
		temp.write("Sweep," + ",".join(to_string(sweep_matrix[run])) + "\n")
		temp.write("\n".join(matrix[run]))
		temp.close()
		os.system("python3 Algorithm_V1.1.py temp_input.csv temp_ouput.csv")
		data_1 = to_float(open("temp_ouput.csv").read().split("\n")[-1].split(",")[1:])
		os.system("rm temp_ouput.csv")
		os.system("rm temp_input.csv")
		data = []
		for i in data_1:
			if i == 0.0:
				data.append(0.000000001)
			else:
	
			data.append(i)		


		## Calculates the divergence of error corrected vs true distribution
		difference = distance.jensenshannon(data,Distribution_in_vitro)

		output_error.append(str(difference))

	outputs.append(output_error)


## Writes out diveregence for the dataset to console. 

for i in outputs:
	print(",".join(i))

for i in Error_range:
	print(str(i))


















