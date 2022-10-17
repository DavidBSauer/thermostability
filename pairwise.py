# taking in an MSA and sequence identity information
#

from matplotlib_venn import venn2
from matplotlib import pyplot as plt
import ETA
import logging
import numpy
import time

def temp_return(txt):
	a = txt.split('|')
	return float(a[-1])


def position_hist_maker(high,low,position,folder):
	alphab = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','-']
	plt.clf()
	plt.close()
	pos = numpy.arange(len(alphab))
	width = 0.4     # gives histogram aspect to the bar diagram
	ax = plt.axes()
	ax.set_xlim(0,21)
	ax.set_xticks(pos + (width / 2))
	ax.set_xticklabels(alphab)
	ax.set_ylim(0,1)
	plt.bar(pos-width/2, high, width, color='r',alpha=0.5,label='Higher Temp Seq')
	plt.bar(pos+width/2, low, width, color='b',alpha=0.5,label='Lower Temp Seq')
#	plt.legend(loc='upper left')
	plt.savefig("./"+folder+"_data/pairwise/figures/histograms/pos_"+str(position+1)+"_histogram.png")
	plt.clf()
	plt.close()


def analyse(MSA,idents,temp,ident,gap,folder):
	used = {}
	length = MSA.get_alignment_length()
	higher_matrix = numpy.zeros(shape=(length,21))
	lower_matrix = numpy.zeros(shape=(length,21))
	ref_matrix = numpy.zeros(shape=(length,1))
	gap_matrix = numpy.zeros(shape=(length,1))
	temp_pairs = 0
	total_pairs = 0
	AA_dict = {'A':0,'C':1,'D':2,'E':3,'F':4,'G':5,'H':6,'I':7,'K':8,'L':9,'M':10,'N':11,'P':12,'Q':13,'R':14,'S':15,'T':16,'V':17,'W':18,'Y':19,'-':20}
	count = 0
	time0 = time.time()	
	for x in MSA:
		count = count+1
		print "pairwise analysis:"+ETA.ETA(count,len(MSA),time0,"1/2N^2")
		used[x.id]=True
		for z in range(0,length):
			if x.seq[z] == "-":
				gap_matrix[z] = gap_matrix[z] +1
		for y in MSA:
			if not(y.id in used):
				if (idents[x.id.split("|")[0]][y.id.split("|")[0]] > ident):
					total_pairs = total_pairs +1
					for z in range(0,length):
						if x.seq[z] != y.seq[z]:
							if (temp_return(x.id)-temp_return(y.id))>temp:
								higher_matrix[z][AA_dict[x.seq[z]]] = higher_matrix[z][AA_dict[x.seq[z]]] +1
								lower_matrix[z][AA_dict[y.seq[z]]] = lower_matrix[z][AA_dict[y.seq[z]]] +1
							elif (temp_return(y.id)-temp_return(x.id))>temp:
								higher_matrix[z][AA_dict[y.seq[z]]] = higher_matrix[z][AA_dict[y.seq[z]]] +1
								lower_matrix[z][AA_dict[x.seq[z]]] = lower_matrix[z][AA_dict[x.seq[z]]] +1

							ref_matrix[z] = ref_matrix[z]+1
					if abs(temp_return(x.id)-temp_return(y.id))>temp:
						temp_pairs = temp_pairs +1
	#calculate 
	logging.info("pairwise: found "+str(temp_pairs)+" pairs that meet the temp condition, of "+str(total_pairs)+" total pairs")
	if (temp_pairs > 0):
		temp_matrix = numpy.sum(higher_matrix,axis=1)/temp_pairs
		ref_matrix = numpy.sum(ref_matrix,axis=1)/total_pairs	
		gap_matrix = numpy.sum(gap_matrix,axis=1)*100/len(MSA)
		higher_matrix = higher_matrix/temp_pairs
		lower_matrix = lower_matrix/temp_pairs	
		diff_matrix = numpy.subtract(temp_matrix,ref_matrix)
				
		numpy.savetxt("./"+folder+"_data/pairwise/higher_residues.csv",higher_matrix,delimiter="\t")
		numpy.savetxt("./"+folder+"_data/pairwise/lower_residues.csv",lower_matrix,delimiter="\t")
		numpy.savetxt("./"+folder+"_data/pairwise/temp_dependent_differences.csv",temp_matrix,delimiter="\t")
		numpy.savetxt("./"+folder+"_data/pairwise/reference_differences.csv",ref_matrix,delimiter="\t")
		numpy.savetxt("./"+folder+"_data/pairwise/temp_dependent_minus_ref_differences.csv",diff_matrix,delimiter="\t")
		numpy.savetxt("./"+folder+"_data/pairwise/gap_frequency.csv",gap_matrix,delimiter="\t")
		for z in range(0,length):
			if (gap_matrix[z] > gap):
				gap_matrix[z] = 0
			else:
				gap_matrix[z] = 1
		diff_matrix_masked = numpy.multiply(diff_matrix,gap_matrix)
		numpy.savetxt("./"+folder+"_data/pairwise/temp_dependent_minus_ref_differences_masked.csv",diff_matrix_masked,delimiter="\t")
		
		plt.clf()

		for x in range (0,length): 
			position_hist_maker(higher_matrix[x],lower_matrix[x],x,folder)
	
		plt.clf()
		width = 1     # gives histogram aspect to the bar diagram
		plt.figure(figsize=(33,18),dpi=300)
		pos = numpy.arange(length)
		ax = plt.axes()
		ax.set_xlim(0,length)
		ax.set_xticks(range(0,length,int(length/10)))
		ax.set_ylim(0,1)
		font_size = 30
		ax.tick_params(axis='both',labelsize=font_size)
		plt.bar(pos,temp_matrix, width, color='r')
		plt.savefig("./"+folder+"_data/pairwise/figures/temp_dependent_differences.png")
		plt.clf()
		plt.close()
		
		width = 1    # gives histogram aspect to the bar diagram
		plt.figure(figsize=(33,18),dpi=300)
		pos = numpy.arange(length)
		ax = plt.axes()
		ax.set_xlim(0,length)
		ax.set_xticks(range(0,length,int(length/10)))
		ax.set_ylim(0,1)
		font_size = 30
		ax.tick_params(axis='both',labelsize=font_size)
		plt.bar(pos,ref_matrix, width, color='b')
		plt.savefig("./"+folder+"_data/pairwise/figures/reference_differences.png")
		plt.clf()
		plt.close()
		
		width = 1     # gives histogram aspect to the bar diagram
		plt.figure(figsize=(33,18),dpi=300)
		pos = numpy.arange(length)
		ax = plt.axes()
		ax.set_xlim(0,length)
		ax.set_xticks(range(0,length,int(length/10)))
		ax.set_ylim(0,1)
		font_size = 30
		ax.tick_params(axis='both',labelsize=font_size)
		plt.bar(pos,diff_matrix, width, color='g')
		plt.savefig("./"+folder+"_data/pairwise/figures/temp_dependent_minus_ref_differences.png")
		plt.clf()
		plt.close()
		
		width = 1     # gives histogram aspect to the bar diagram
		plt.figure(figsize=(33,18),dpi=300)
		pos = numpy.arange(length)
		ax = plt.axes()
		ax.set_xlim(0,length)
		ax.set_xticks(range(0,length,int(length/10)))
		ax.set_ylim(0,1)
		font_size = 30
		ax.tick_params(axis='both',labelsize=font_size)
		plt.bar(pos,diff_matrix_masked, width, color='g')
		plt.savefig("./"+folder+"_data/pairwise/figures/temp_dependent_minus_ref_differences_gap_positions_masked.png")
		plt.clf()
		plt.close()

