from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
import logging
from matplotlib_venn import venn2
from matplotlib import pyplot as plt

def spec(txt):
	a = txt.split('|')
	b = a [1]
	return b.lower()
	
def sort(MSA,species_temp,out_file):

	spec2 = []
	temps = []
	spec_used ={}
	for x in MSA:
		if spec(x.id) in species_temp:
			spec2.append(SeqRecord(Seq(str(x.seq),x.seq.alphabet), x.id+"|"+str(species_temp[spec(x.id)]),'',''))
			temps.append(species_temp[spec(x.id)])
			spec_used [spec(x.id)]= True

	MSA2 = MultipleSeqAlignment(spec2)		
	MSA2.sort(key = lambda record: species_temp[spec(record.id)])
	AlignIO.write(MSA2, "./"+out_file+"_data/"+out_file+"_sorted.fa", "fasta")
	logging.info("of "+str(len(MSA))+" sequences, and "+str(len(species_temp.keys()))+" species with temp data "+str(len(MSA2))+" were in both lists and saved")
	plt.close('all')
	v1 = venn2(subsets=(len(MSA),0,len(MSA2)),set_labels = ('Input sequences', '', 'Assigned seqs')) #input seq vs used seq
	v1.get_label_by_id("01").set_text("")
	v1.get_patch_by_id('10').set_color('red')
	v1.get_patch_by_id('11').set_alpha(0.5)
	v1.get_patch_by_id('11').set_color('blue')
	plt.savefig("./"+out_file+"_data/"+out_file+"_sequences_assigned_of_input_venn.png")
	plt.close('all')
	v2 = venn2(subsets=(len(species_temp.keys()),0,len(spec_used)),set_labels = ('Species data', '', 'Used spec data')) #input spec vs used species
	v2.get_label_by_id("01").set_text("")
	v2.get_patch_by_id('10').set_color('green')
	v2.get_patch_by_id('11').set_alpha(0.5)
	v2.get_patch_by_id('11').set_color('blue')
	plt.savefig("./"+out_file+"_data/"+out_file+"_species_temp_data_used_venn.png")	
	plt.close('all')
	plt.figure(figsize=(10,8), dpi=300)
	plt.hist(temps, range(-10,121,5), normed=False, facecolor='green', alpha=0.5)
	plt.tick_params(axis='x', which='both',bottom='on',top='off')
	plt.xlabel('OGT')
	plt.ylabel('Count')
	plt.tight_layout()
	plt.savefig("./"+out_file+"_data/temps_assigned.png")
	plt.close('all')

	return MSA2

def assign(MSA,species_temp,out_file):

	spec2 = []
	temps = []
	spec_used ={}
	for x in MSA:
		if spec(x.id) in species_temp:
			spec2.append(SeqRecord(Seq(str(x.seq),x.seq.alphabet), x.id,'',''))
			temps.append(species_temp[spec(x.id)])
			spec_used [spec(x.id)]= True
	MSA2 = MultipleSeqAlignment(spec2)		
	MSA2.sort(key = lambda record: species_temp[spec(record.id)])
	AlignIO.write(MSA2, "./"+out_file+"_data/"+out_file+"_sorted.fa", "fasta")
	logging.info("of "+str(len(MSA))+" sequences, and "+str(len(species_temp.keys()))+" species with temp data "+str(len(spec2))+" were in both lists and saved")
	plt.close('all')
	v1 = venn2(subsets=(len(MSA),0,len(spec2)),set_labels = ('Input Seq.', '', 'Assigned seqs')) #input seq vs used seq
	v1.get_label_by_id("01").set_text("")
	v1.get_patch_by_id('10').set_color('red')
	v1.get_patch_by_id('11').set_alpha(0.5)
	v1.get_patch_by_id('11').set_color('blue')
	plt.savefig("./"+out_file+"_data/"+out_file+"_sequences_assigned_of_input_venn.png")
	plt.close('all')
	v2 = venn2(subsets=(len(species_temp.keys()),0,len(spec_used)),set_labels = ('Species data', '', 'Used species data')) #input spec vs used species
	v2.get_label_by_id("01").set_text("")
	v2.get_patch_by_id('10').set_color('green')
	v2.get_patch_by_id('11').set_alpha(0.5)
	v2.get_patch_by_id('11').set_color('blue')
	plt.savefig("./"+out_file+"_data/"+out_file+"_species_temp_data_used_venn.png")
	plt.close('all')
	plt.figure(figsize=(10,8), dpi=300)
	plt.hist(temps, range(-10,121,5), normed=False, facecolor='green', alpha=0.5)
	plt.tick_params(axis='x', which='both',bottom='on',top='off')
	plt.xlabel('OGT')
	plt.ylabel('Count')
	plt.tight_layout()
	plt.savefig("./"+out_file+"_data/temps_assigned.png")
	plt.close('all')
