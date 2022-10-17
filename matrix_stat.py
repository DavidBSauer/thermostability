# script to take in an MSA where the header includes the organism growth temperature
# split into clades based on temperature
# then calculate the frequency of amino acid by position for each clade
# and subtract these frequencies 
# print out csv files for each set and 

import numpy
from matplotlib import cm
import matplotlib.pyplot as plt
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
import logging

 

def temp_return(txt):
	a = txt.split('|')
	return a[-1]

def figure(matrix,title,folder):
	plt.clf()
	plt.close()
	new_matrix = matrix.transpose()
	new_matrix = numpy.flipud(new_matrix)
	(x,length) = new_matrix.shape
	plt.figure(figsize=(33,18), dpi=300)
	ax = plt.axes()
	plt.tick_params(axis='y', which='both',left='off',right='off')
	yticklabels=['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','-']
	pos = numpy.arange(len(yticklabels))
	pos = pos[::-1]
	ax.set_yticks(pos + (1.0 / 2))
	xticks = range(0,length+1,100)
	xticks.pop(0)
	ax.set_xticks(xticks)
	ax.set_yticklabels(yticklabels)
	font_size = 30
	ax.tick_params(axis='both',labelsize=font_size)
	if (numpy.amin(matrix)>=0):
		p = ax.pcolormesh(new_matrix)
	else:
		p = ax.pcolormesh(new_matrix,cmap=cm.coolwarm)
	cbar = plt.colorbar(p)
	cbar.ax.tick_params(labelsize = font_size)
	plt.tight_layout()
	plt.savefig("./"+folder+"_data/global/figures/"+title+'.png')
	plt.clf()
	plt.close()

def hist_maker(temps,folder):
	plt.clf()
	plt.close()
	plt.figure(figsize=(10,8), dpi=300)
	plt.hist(temps, range(-10,121,5), normed=False, facecolor='green', alpha=0.5)
	plt.tick_params(axis='x', which='both',bottom='on',top='off')
	plt.xlabel('OGT')
	plt.ylabel('Count')
	plt.tight_layout()
	plt.savefig("./"+folder+"_data/global/figures/global_temps_used_histogram.png")
	plt.clf()
	plt.close()

def position_hist_maker(high,low,position,folder):
	alphab = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','-']
	plt.clf()
	plt.close()
	plt.figure(figsize=(10,8), dpi=300)
	pos = numpy.arange(len(alphab))
	width = 0.4    # gives histogram aspect to the bar diagram
	ax = plt.axes()
	ax.set_xlim(-width,21)
	ax.set_xticks(pos + (width / 2))
	ax.set_xticklabels(alphab)
	ax.set_ylim(0,1)
	plt.tick_params(axis='x', which='both',bottom='on',top='off')
	plt.bar(pos-width/2, high, width, color='r',alpha=0.5,label='Thermo+Hyperthermophiles')
	plt.bar(pos+width/2, low, width, color='b',alpha=0.5,label='Mesophiles')
	plt.tight_layout()
	plt.savefig("./"+folder+"_data/global/figures/histograms/pos_"+str(position+1)+"_histogram.png")
	plt.clf()
	plt.close()


def divide(MSA,folder):
	length = MSA.get_alignment_length()

	thermo_matrix = numpy.zeros(shape=(length,21))
	hyper_matrix = numpy.zeros(shape=(length,21))
	meso_matrix = numpy.zeros(shape=(length,21))
	psych_matrix = numpy.zeros(shape=(length,21))

	# make dictionary of AA to number
	AA_dict = {'A':0,'C':1,'D':2,'E':3,'F':4,'G':5,'H':6,'I':7,'K':8,'L':9,'M':10,'N':11,'P':12,'Q':13,'R':14,'S':15,'T':16,'V':17,'W':18,'Y':19,'-':20}

	#go through each matrix, take a slice and count all AAs
	temps =[]
	thermo_entries = []
	hyper_entries  = []
	meso_entries = []
	psychro_entries = []
	

	for x in MSA:
		seq_temp = float(temp_return(x.id))
		temps.append(seq_temp)
		if (seq_temp >= 60):
			for y in range(0, length):
				hyper_matrix [y] [AA_dict[x.seq[y]]] += 1
			hyper_entries.append(SeqRecord(Seq(str(x.seq),x.seq.alphabet), x.id,'',''))
		elif ((seq_temp <60) & (seq_temp >= 45)):
			for y in range(0, length):
				thermo_matrix [y] [AA_dict[x.seq[y]]] += 1
			thermo_entries.append(SeqRecord(Seq(str(x.seq),x.seq.alphabet), x.id,'',''))
		elif ((seq_temp <45) & (seq_temp >=20)):
			for y in range(0, length):
				meso_matrix [y] [AA_dict[x.seq[y]]] += 1
			meso_entries.append(SeqRecord(Seq(str(x.seq),x.seq.alphabet), x.id,'',''))
		else: 
			for y in range(0, length):
				psych_matrix [y] [AA_dict[x.seq[y]]] += 1
			psychro_entries.append(SeqRecord(Seq(str(x.seq),x.seq.alphabet), x.id,'',''))

	hyper_thermo = (hyper_matrix + thermo_matrix)

	#normalize the matrices to the number of sequences, get frequencies
	if (len(hyper_entries)>0):
		hyper_matrix_norm = hyper_matrix / len(hyper_entries)
		numpy.savetxt("./"+folder+"_data/global/hyper_norm.csv",hyper_matrix_norm,delimiter="\t")
		numpy.savetxt("./"+folder+"_data/global/hyper.csv",hyper_matrix,delimiter="\t")
		figure(hyper_matrix_norm,'hyper_norm',folder)
		hyper_MSA = MultipleSeqAlignment(hyper_entries)
		AlignIO.write(hyper_MSA, "./"+folder+"_data/global/hyperthermophiles.fa", "fasta")
		logging.info("the number of hyperthermophiles: "+str(len(hyper_entries)))
	if (len(thermo_entries)>0):
		thermo_matrix_norm = thermo_matrix / len(thermo_entries)
		numpy.savetxt("./"+folder+"_data/global/thermo_norm.csv",thermo_matrix_norm,delimiter="\t")
		numpy.savetxt("./"+folder+"_data/global/thermo.csv",thermo_matrix,delimiter="\t")
		figure(thermo_matrix_norm,'thermo_norm',folder)
		thermo_MSA = MultipleSeqAlignment(thermo_entries)
		AlignIO.write(thermo_MSA, "./"+folder+"_data/global/thermophiles.fa", "fasta")
		logging.info("the number of thermophiles: "+str(len(thermo_entries)))
	if (len(meso_entries) >0):
		meso_matrix_norm = meso_matrix / len(meso_entries)
		numpy.savetxt("./"+folder+"_data/global/meso_norm.csv",meso_matrix_norm,delimiter="\t")
		numpy.savetxt("./"+folder+"_data/global/meso.csv",meso_matrix,delimiter="\t")
		figure(meso_matrix_norm,'meso_norm',folder)
		meso_MSA = MultipleSeqAlignment(meso_entries)
		AlignIO.write(meso_MSA, "./"+folder+"_data/global/mesophiles.fa", "fasta")
		logging.info("the number of mesophiles: "+str(len(meso_entries)))
	if (len(psychro_entries)>0):
		psych_matrix_norm = psych_matrix / len(psychro_entries)
		numpy.savetxt("./"+folder+"_data/global/psych_norm.csv",psych_matrix_norm,delimiter="\t")
		numpy.savetxt("./"+folder+"_data/global/psych.csv",psych_matrix,delimiter="\t")
		figure(psych_matrix_norm,'psych_norm',folder)
		psychro_MSA = MultipleSeqAlignment(psychro_entries)
		AlignIO.write(psychro_MSA, "./"+folder+"_data/global/psychrophiles.fa", "fasta")
		logging.info("the number of psychrophiles: "+str(len(psychro_entries)))
	if ((len(hyper_entries) + len(thermo_entries))>0):
		hyper_thermo = hyper_matrix + thermo_matrix	
		hyper_thermo_norm = (hyper_matrix + thermo_matrix) / (len(hyper_entries) + len(thermo_entries))		
		numpy.savetxt("./"+folder+"_data/global/hyper+thermo_norm.csv",hyper_thermo_norm,delimiter="\t")
		figure(hyper_thermo_norm,'hyper_thermo_norm',folder)

	hist_maker(temps,folder)
	# take differences
	if ((len(hyper_entries)>0) and (len(meso_entries)>0)):
		hyper_meso_norm = hyper_matrix_norm - meso_matrix_norm
		numpy.savetxt("./"+folder+"_data/global/hyper-meso_norm.csv",hyper_meso_norm,delimiter="\t")
		figure(hyper_meso_norm,'hyper_meso_diff',folder)
	if (((len(hyper_entries) +len(thermo_entries)) >0) and (len(meso_entries)>0)):
		hyper_thermo_meso = hyper_thermo_norm - meso_matrix_norm
		numpy.savetxt("./"+folder+"_data/global/hyper+thermo-meso_norm.csv",hyper_thermo_meso,delimiter="\t")
		figure(hyper_thermo_meso,'hyper_thermo_meso_diff',folder)
		p_value = numpy.ones(shape=(length,21))
		for x in range (0,length): 
			position_hist_maker(hyper_thermo_norm[x],meso_matrix_norm[x],x,folder)

