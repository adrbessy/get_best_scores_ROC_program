### ----------------------------
### Best scores ROC program
### ----------------------------

'''
This program allows to calculate the area under the ROC (receiver operating characteristic) curve between bound sequences and unbound sequences.
It takes into account the best score of each sequence.
It can compare several matrices.
So you need _a matrix or several matrix with frequency values , _fasta bound sequences, __fasta unbound sequences
This program was written by Adrien Bessy and Arnaud Stigliani with the collaboration of Francois Parcy and was inspired by Morpheus program written by Eugenio Gomez Minguet
and (4-scores.py and CalculateScoreAboveTh_and_Interdistances.py) programs written by Laura Gregoire.
'''

import numpy as np
from Bio import SeqIO
import time
import sys
#from tqdm import *  
from math import exp
from math import log
from datetime import datetime
import decimal
from decimal import Decimal    
import argparse
import re

parser = argparse.ArgumentParser()                                               

parser.add_argument("--factor", "-fac", type = str, default= "ARF2")
parser.add_argument("--pseudoCount", "-pc",type = float, default = 0.001)
parser.add_argument("--sequence_number", "-sequence_number", type = int, default= 11654)

args = parser.parse_args()

#python get_best_scores_ROC.py -fac "ARF2" -pc 0.001 -sequence_number 11654
#python get_best_scores_ROC.py -fac "ARF5" -pc 0.001 -sequence_number 26659


factorTranscription = args.factor
pseudoCount = args.pseudoCount
sequence_number = args.sequence_number

if factorTranscription == "ARF2" :
	FastaFile = "../sequences/ARF2_bound_sequences.fas"  
	MatrixFile1 = "../matrices/ARF2_OMalley_matrixC.txt" 
	matrixType1 = "freq" 
	dependencyFile = ""
	MatrixFile2 = "../get_base_frequencies/ARF2_1500Seq_5prime_freq_pasteTo'OMalley.txt" 
	matrixType2 = "freq" 
	FastaFileN = "../sequences/ARF2_neg1.fas"
	dependencyFile2 = ""
	#MatrixFile3 = "../get_base_frequencies/ARF2_1500_5prime_freq.txt" 
	#matrixType3 = "freq" 
	#dependencyFile3 = ""
	
	
if factorTranscription == "ARF5" :
	FastaFile = "../sequences/ARF5_bound_sequences.fas"  
	MatrixFile1 = "../matrices/ARF5_OMalley_matrixC.txt" 
	matrixType1 = "freq" 
	dependencyFile = ""
	MatrixFile2 = "../get_base_frequencies/ARF5_1500Seq_3prime3bases_freq_pasteTo'OMalley.txt" 
	matrixType2 = "freq" 
	dependencyFile2 = ""
	FastaFileN = "../sequences/ARF5_neg1.fas"
	#MatrixFile3 = "../get_base_frequencies/ARF5_1500_3prime_freq.txt" 
	#matrixType3 = "freq" 
	#dependencyFile3 = ""


''' separation between numbers can be spaces, tabulation, comas...
                                                                         Example :   A C G T
                                                                  position 1:           0.16456   0.21614       0.1565,0.1645
                                                                  position 2:           0.645; 0.654    0.155 |||||| 0.4444
                                                                                        ...
'''
codigo = { 'AA': 0, 'AC': 1, 'AG': 2, 'AT': 3,
		'CA': 4, 'CC': 5, 'CG': 6, 'CT': 7,
		'GA': 8, 'GC': 9, 'GG': 10, 'GT': 11,
		'TA': 12, 'TC': 13, 'TG': 14, 'TT': 15 }
           
codigoi = { "A" : "T", "C" : "G", "G" : "C", "T" : "A"}


def seq_c(site):
        site_i = site[-1::-1]
        site_c = ""
        for x in site_i:
		y = codigoi[x]
                site_c = site_c + y
	return site_c 

def get_dependency_matrix(dependencyFile,num) : 
	G = open(dependencyFile,"r")
	dependency_file_content = G.read().replace("\r","\n") + "\n"
	G.close()
	
	num2 = re.compile(r"(?<![\d.])[0-9]+(?![\d.])")
	position_dependency_matrix = num2.findall(dependency_file_content)
	position_dependency_matrix = map(int, position_dependency_matrix)
	
	dependency_matrix = num.findall(dependency_file_content)
	dependency_matrix = map(float, dependency_matrix)
	
	dependency_matrix_associated_with_position = []
	index1 = 0
	index2 = 3
	index3 = 0
	index4 = 64
	for i in range(0, 3):
		dependency_matrix_associated_with_position.append(position_dependency_matrix[index1:index2])
		dependency_matrix_associated_with_position.append(dependency_matrix[index3:index4])
		index1 = index1 + 3
		index2 = index2 + 3
		index3 = index3 + 64
		index4 = index4 + 64
		
	return(dependency_matrix_associated_with_position)

def add_scores_associated_with_interdependent_positions(dependency_matrix,scoreStrandPos,scoreStrandNeg,strandPos):
	cStrand = ""
	for lettre in seq_c(strandPos):
		cStrand = lettre + cStrand
	cStrand = cStrand[::-1]
	
	site1 = ""
	Csite1 = ""
	for i in dependency_matrix[0]:
		site1 = site1 + strandPos[i-1]
		Csite1 = Csite1 + cStrand[i-1]
	scoreStrandPos = scoreStrandPos + dependency_matrix[1][codigo[site1]]
	scoreStrandNeg = scoreStrandNeg + dependency_matrix[1][codigo[Csite1]]
	return(scoreStrandPos, scoreStrandNeg)

def calculateBestScores(MatrixFile1,FastaFile,matrixType,dependencyFile) :        

	###first matrix###
	# These 3 lines allows to retrieve the matrix from the file
	F = open(MatrixFile1,"r")
	matrix = F.read().replace("\r","\n") + "\n"
	F.close()

	# These 3 lines allows to retrieve all the individual frequency values from the matrix and put them in order into a list
	import re
	num = re.compile(r"([+-]?\d+[.,]\d+)")
	regex = re.compile('[0-9]+')
	

	## These lines allows to transform the frequency values into scores values
	freq_mat = []
	if matrixType == "count" :
		Mdata = regex.findall(matrix)
		for i in range(0,len(Mdata)):
			if i%4==0:
				#print("Mdata[i] : ",Mdata[i])
				Sum = float(Mdata[i]) + float(Mdata[i+1]) + float(Mdata[i+2]) + float(Mdata[i+3])
				#print("Sum : ",Sum)
				one = float(Mdata[i]) / Sum
				two = float(Mdata[i+1]) / Sum
				three = float(Mdata[i+2]) / Sum
				four = float(Mdata[i+3]) / Sum
				freq_mat.append(one)
				freq_mat.append(two)
				freq_mat.append(three)
				freq_mat.append(four)

		matF = []
		lenMotif=0
		for i in range(0,len(freq_mat)):
			if i%4==0:
				lenMotif=lenMotif+1
				fmax = float(max(freq_mat[i],freq_mat[i+1],freq_mat[i+2],freq_mat[i+3])) + pseudoCount
				for j in range (0,4):
					matF.append(np.log(float(float(freq_mat[i+j]) + pseudoCount) /fmax))
	if matrixType == "freq" :
		Mdata = num.findall(matrix)
		matF = []
		lenMotif=0
		for i in range(0,len(Mdata)):
			if i%4==0:
				lenMotif=lenMotif+1
				fmax = float(max(Mdata[i],Mdata[i+1],Mdata[i+2],Mdata[i+3])) + pseudoCount
				for j in range (0,4):
					matF.append(np.log(float(float(Mdata[i+j]) + pseudoCount) /fmax))
			
	#print("matF : ",matF)			
	matRev = list(reversed(matF))
	
	# This line allows to retrieve all the sequences from the fasta file
	print("FastaFile : ",FastaFile) 
	sequences = SeqIO.to_dict(SeqIO.parse(FastaFile, "fasta"))

	print "  Il y a %s sequence(s) a analyser"%(sequence_number)

	# The following line allows to produce the reversed matrix
	'''if we take the example given before : A T G C
				Position 1:      0.4444  0.155  0.654   0.645
				Position 2:      0.1645  0.1565 0.21614 0.16456
	Now, we can notice that scores change between the positions 1 and 2, and between A and T, and between G and C.
	So we can calculate with this reverse matrix, the score of the complementary strand.
	'''

	# We will store in this list all the interdistances between motifs found in all sequences.

	# we write a file where will be stored every scores
	fileScores = open("%s.scores"%FastaFile,"w")
	fileScores.write("name\tposition\tsens\tscore\n")
	i=1
	 
	# here we store the best score for each sequence
	bestScoreBySeq=[]
	nb = 0
	# We look at all the fasta sequences:
	#for s in tqdm(sequences):
	for s in sequences:
		# We will store in this list all the best scores (see the threshold after) found for subsequences of one sequence
		GoodScorePositionsStrandPos=[]

		# This line allows to retrieve the DNA sequence
		seq=sequences[s].seq
		score_seq = []
		# We look at each sub-sequences of the whole sequence. Each sub-sequence has the same length that the matrix length.
		for c in range(len(seq)-lenMotif+1):
			strandPos = seq[c:c+lenMotif].upper()
			test = 0
			for nu in strandPos :
				if nu not in ["A","C","G","T"]:
					test=1
			if test == 1:
				score = "NA"
			else :
				n=0
				#These lines allows to calculate a score for one sub-sequence
				scoreStrandPos=0
				scoreStrandNeg=0
				while n<lenMotif:
					if strandPos[n] == 'A':
						scoreStrandPos=scoreStrandPos +matF[n*4]
						scoreStrandNeg=scoreStrandNeg +matRev[n*4]
					elif strandPos[n] == 'C':
						scoreStrandPos=scoreStrandPos +matF[n*4+1]
						scoreStrandNeg=scoreStrandNeg +matRev[n*4+1]
					elif strandPos[n] == 'G':
						scoreStrandPos=scoreStrandPos +matF[n*4+2]
						scoreStrandNeg=scoreStrandNeg +matRev[n*4+2]
					elif strandPos[n] == 'T':
						scoreStrandPos=scoreStrandPos +matF[n*4+3]
						scoreStrandNeg=scoreStrandNeg +matRev[n*4+3]  
					n += 1
				if dependencyFile != "" : 			
					scoreStrandPos, scoreStrandNeg = add_scores_associated_with_interdependent_positions(get_dependency_matrix(dependencyFile,num),scoreStrandPos,scoreStrandNeg,strandPos)
				# We just keep the maximum score value between the plus strand and the negative strand
				scores = max (scoreStrandPos, scoreStrandNeg)
				#print("scores : ",scores)
				score_seq.append(scores)
				#print("score_seq : ",score_seq.append(scores))
		bestScoreBySeq.append(max(score_seq))
		if sequence_number :
			nb = nb + 1
		if nb == sequence_number : 
			break
			
	return bestScoreBySeq
                                    
def ROC(data1, data2, Min, Max):
	output = MatrixFile1+"Best-scores_ROC-AUC.txt"
	R = open(output, "w")
	R.write("@Rango\tFreqData1(y)\tFreqData2(x)\tCountData1\tCountData2\n")
	l1 = len(data1)
	l2 = len(data2)
	A = (Max-Min) / 10000
	limit = Max
	Freqx = []
	Freqy = []
	while limit > Min:
		c1 = 0
		limit = limit - A
		for i in data1:
			if i < limit:
				continue
			else:
				c1 = c1 + 1
		c2 = 0
		for h in data2:
			if h < limit:
				continue
			else:
				c2 = c2 + 1
		prop1 = Decimal(c1) / Decimal(l1)
		prop2 = Decimal(c2) / Decimal(l2)
		f1 = prop1.quantize(Decimal("1.00000000"))
		f2 = prop2.quantize(Decimal("1.00000000"))
		R.write(str(limit) + "\t" + str(f1) + "\t" + str(f2) + "\t" + str(c1) + "\t" + str(c2) + "\n")
		Freqx.append(prop2)
		Freqy.append(prop1)
	## R calculation
	##   x value -> frequency seqs NotBound  (Freqx)
	##   y value -> frequency seqs Bound     (Freqy)
	## Trapeze method
	i = 0
	AreaT = Decimal(0)
	l = len(Freqx) - 1
	while i < l:
		I = Freqx[i+1] - Freqx[i]
		if Freqy[i] == 0:
			i = i + 1
			continue
		if I == 0:
			i = i + 1
			continue
		else:
			A = I * (Freqy[i] + Freqy[i+1])/2
			AreaT = AreaT + A
			i = i + 1
	AreaT = str(round(float(AreaT),5))
	R.write("\nR (ROC-AUC) = " + AreaT + "\n")
	R.close()
	return AreaT,output,Freqx,Freqy
              
data1 = calculateBestScores(MatrixFile1,FastaFile,matrixType1,dependencyFile)
data2 = calculateBestScores(MatrixFile1,FastaFileN,matrixType1,dependencyFile)
data3 = calculateBestScores(MatrixFile2,FastaFile,matrixType2,dependencyFile2)
data4 = calculateBestScores(MatrixFile2,FastaFileN,matrixType2,dependencyFile2)
#data5 = calculateBestScores(MatrixFile3,FastaFile,matrixType3,dependencyFile3)
#data6 = calculateBestScores(MatrixFile3,FastaFileN,matrixType3,dependencyFile3)

Max1 = max(max(data1),max(data2))
Min1 = min(min(data1),min(data2))
Max2 = max(max(data3),max(data4))
Min2 = min(min(data3),min(data4))
#Max3 = max(max(data5),max(data6))
#Min3 = min(min(data5),min(data6))
    
AreaT1,output,Freqx1,Freqy1 = ROC(data1, data2, Min1, Max1)
AreaT2,output,Freqx2,Freqy2 = ROC(data3, data4, Min2, Max2)
#AreaT3,output,Freqx3,Freqy3 = ROC(data5, data6, Min3, Max3)

############################
## Graphic representation ##
############################

import numpy as np
import matplotlib.pyplot as plot

print "\nGenerating image files..."

ROCvalue1 = MatrixFile1+" ROC-AUC = " + AreaT1
ROCvalue2 = MatrixFile2+" ROC-AUC = " + AreaT2
#ROCvalue3 = MatrixFile3+" ROC-AUC = " + AreaT3

plot.title('ROC-AUC of the best scores of each sequence with PWMs')
plot.plot(Freqx1, Freqy1, marker="None", linestyle="-", color="purple", lw="3")
plot.plot(Freqx2, Freqy2, marker="None", linestyle="-", color="cornflowerblue", lw="3")
#plot.plot(Freqx3, Freqy3, marker="None", linestyle="-", color="brown", lw="3")
plot.axis([0,1,0,1])
plot.xlabel("unbound sequences")
plot.ylabel("bound sequences")
plot.figtext(0.45,0.45, ROCvalue1,color="purple")
plot.figtext(0.45,0.35, ROCvalue2,color="cornflowerblue")
#plot.figtext(0.45,0.25, ROCvalue3,color="brown")

plot.show()              
