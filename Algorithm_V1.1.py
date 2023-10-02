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
# 			python3 Input.csv Output.csv
#
# 
# Changing various script inputs, outputs and parameters can be done as follows
#
# Tab seperate input and output: replace all .split(",") with .split("\t") 
# Change number of iterations: replace iterations = 100 with iterations = required number
# Take sweep data as a fraction : replace Sweep.append(i/100) with Sweep.append(i)
# Include error correction : Uncomment in def run_algorithm function.
#
# ------------------------------------------------------------------------------------------------


import os
import sys
import itertools
import random
import math



def phasome(A):
	E = []
	for B in A:
		C = ""
		for D in B:
			C = C + str(D)
		E.append(C)
	return E

def normalise(A):
	B = []
	if sum(A) != 0:
		for i in A:
			B.append(i/sum(A))
	else:
		B = A
	return B

def to_string(A):
	B = []
	for i in A:
		B.append(str(i))
	return B

def to_float(A):
	B = []
	for i in A:
		B.append(float(i))
	return B

def ON_Ratio(dist,phasomes):
	ON_Dist = []

	for gene in range(0,len(phasomes[0])):
		pos = 0			
		ON = 0
		OFF = 0
		for p in phasomes:
			if p[gene] == "1":
				ON = ON + dist[pos]
			pos = pos + 1
		ON_Dist.append(ON)

	return ON_Dist


def manhattan(a, b):
	var = 0
	for i in range(0,len(a)):
		if a[i] == b[i]:
			var = var + 1
	return var


def gene_adjusted(previous_sweep,Sweep_Aim):
	squared_diff = []
	pos = 0 
	
	for gene in previous_sweep:
		squared_diff.append(math.sqrt((gene - Sweep_Aim[pos])**2))

		pos = pos + 1

	adjustment =  Sweep_Aim[squared_diff.index(max(squared_diff))] - previous_sweep[squared_diff.index(max(squared_diff))]

	direction = ""
	if adjustment < 0:
		direction = "decrease"
	else:
		direction = "increase"

	output = [squared_diff.index(max(squared_diff)), adjustment, direction]

	return output






def genes_to_swap(corrected_distributions, Phasotypes, correction):

	current_phasotypes = []

	pos = 0
	for i in corrected_distributions:
		if 0 < i:
			current_phasotypes.append(Phasotypes[pos])
		pos = pos + 1


	present_adjust = []
	to_adjust = []

	for i in current_phasotypes:
		if correction[-1] == "increase":
			if i[correction[0]] == "0":
				present_adjust.append(i)
				to_adjust.append(i[:correction[0]] + "1" + i[correction[0] + 1:])	
		if correction[-1] == "decrease":
			if i[correction[0]] == "1":
				present_adjust.append(i)
				to_adjust.append(i[:correction[0]] + "0" + i[correction[0] + 1:])

	return present_adjust, to_adjust

def normalise_adjustment(present_adjust,to_adjust,previous_distribution,Phasotypes,correction):

	value_present = []
	pos = 0
	for p in Phasotypes:
		if p in present_adjust:
			value_present.append(previous_distribution[pos])

		pos = pos + 1

	normalised_value_present = []

	for i in normalise(value_present):
		normalised_value_present.append(correction[1]*i)
	
	return normalised_value_present

def correct_dist(present_adjust,to_adjust,previous_distribution,Phasotypes,correction,normalised_value_present):

	adjusted_distribution = []
	
	pos = 0 

	for phase in Phasotypes:
		if phase in present_adjust:
			if correction[-1] == "increase":
				new_value = previous_distribution[Phasotypes.index(phase)] - normalised_value_present[present_adjust.index(phase)]
				adjusted_distribution.append(new_value)
			elif correction[-1] == "decrease":
				new_value = previous_distribution[Phasotypes.index(phase)] + normalised_value_present[present_adjust.index(phase)]
				adjusted_distribution.append(new_value)

		elif phase in to_adjust:
			if correction[-1] == "increase":
				new_value = previous_distribution[Phasotypes.index(phase)] + normalised_value_present[to_adjust.index(phase)]
				adjusted_distribution.append(new_value)
			elif correction[-1] == "decrease":
				new_value = previous_distribution[Phasotypes.index(phase)] - normalised_value_present[to_adjust.index(phase)]
				adjusted_distribution.append(new_value)

		else:
			adjusted_distribution.append(previous_distribution[Phasotypes.index(phase)])

	return adjusted_distribution




def parsimony_algorithm(Distribution, Sweep_Aim, Phasotypes):

	corrected_sweeps = [ON_Ratio(Distribution,Phasotypes)]
	corrected_distributions = [Distribution]
	iterations = 100


	for iteration in range(0, iterations):
		previous_sweep = corrected_sweeps[-1]
		previous_distribution = corrected_distributions[-1]

		correction = gene_adjusted(previous_sweep,Sweep_Aim)

		present_adjust, to_adjust = genes_to_swap(previous_distribution, Phasotypes, correction)
	
		normalised_value_present = normalise_adjustment(present_adjust,to_adjust,previous_distribution,Phasotypes,correction)
		
		adjusted_distribution = correct_dist(present_adjust,to_adjust,previous_distribution,Phasotypes,correction,normalised_value_present)


		corrected_distributions.append(adjusted_distribution)
		corrected_sweeps.append(ON_Ratio(adjusted_distribution,Phasotypes))

		if adjusted_distribution == corrected_distributions[-1] and adjusted_distribution == corrected_distributions[-2]:
			break

		for i in adjusted_distribution:
			if i < -0.01 or 1.01 < i:
				print("ERROR outside of frame detected in iteration: " + str(iteration))
				quit()
				break
		if sum(adjusted_distribution) < -0.01 or 1.01 < sum(adjusted_distribution):
			print("ERROR outside of frame detected in iteration: " + str(iteration))
			quit()



	return corrected_sweeps, corrected_distributions



def Find_Missing_Data(Dist_Input, Missing_Data, Phasotypes):

	threshold = 200
	
	threshold_missing = []

	for i in Missing_Data:
		if len(i) != 0 and i.count("N") < len(i)*threshold:
			threshold_missing.append(i)

	Ratios = []
	for i in normalise(Dist_Input):
		Ratios.append(i + 1)
	corrected_missing = []

	for value in threshold_missing:
		Distances = []
		for p in Phasotypes:
			Distances.append(manhattan(value,p))

		weighted = []
		
		for d in range(0,len(Distances)):
			weighted.append(Distances[d]*Ratios[d])

		corrected_missing.append(Phasotypes[weighted.index(max(weighted))])

	for correct in corrected_missing:
		Dist_Input[Phasotypes.index(correct)] = Dist_Input[Phasotypes.index(correct)] + 1

	return Dist_Input


def run_algorithm(file, output_file):


	Input_file = open(file).read().split("\n")

	Phasotype = phasome(list(itertools.product([0, 1], repeat=len(Input_file[0].split(",")[1:]))))

	Sweep = []
	for i in to_float(Input_file[1].split(",")[1:]):
		if i == "0":
			Sweep.append(1/100000000)
		else:
			Sweep.append(i/100)


	Distribution = [0]*len(Phasotype)
	Missing  = []
	for line in Input_file[2:]:
		l = line.split(",")
		phasotype = "".join(l[1:])
		if phasotype in Phasotype:
			Distribution[Phasotype.index(phasotype)] = Distribution[Phasotype.index(phasotype)] + 1
		elif phasotype not in Phasotype and len(phasotype) != 0:
			Missing.append(phasotype)



	## This is an optional call -- if you want to include and correct for missing data using parsimony un-comment (remove the """ from the top and bottom)


	if len(Missing) != 0:
		error_correct_distribution = Find_Missing_Data(Distribution, Missing,Phasotype)
		Distribution = error_correct_distribution


	ratio_of_phasotypes = normalise(Distribution)
	corrected_sweeps, corrected_distributions = parsimony_algorithm(ratio_of_phasotypes,Sweep, Phasotype)

	Output_file = open(output_file,"w")
	Output_file.write("Corrected Colony Sweep," + ",".join(to_string(corrected_sweeps[-1])) + "\n") 
	Output_file.write("Phasotypes," + ",".join(Phasotype) + "\n")
	Output_file.write("Corrected Distribution," + ",".join(to_string(corrected_distributions[-1])))
	Output_file.close()



if len(sys.argv) == 3:
	run_algorithm(sys.argv[1], sys.argv[2])
else:
	print("Please provide a valid Input file and Output File")























