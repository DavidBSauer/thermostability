from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import logging
import random
from Bio import AlignIO
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from matplotlib_venn import venn2
from matplotlib import pyplot as plt

def header(txt):
	return txt.split('|')[0]

def cleanup(MSA,new_file):

	seq_num = len(MSA)
	length = len(MSA[-1])

	seqs =[]
	for x in MSA:
		if ('X' not in str(x.seq).upper()):
			new_seq = str(x.seq).upper()
			seqs.append(SeqRecord(Seq(new_seq,x.seq.alphabet), x.id,'',''))

	seq_num2 = len(seqs)
	SeqIO.write(seqs, "./"+new_file+"_data/"+new_file+"_Xed.fa", "fasta")
	logging.info(str("of " + str(seq_num) + " sequences starting, " + str(seq_num2) + " sequences were kept"))
	return seqs



def cleanup_already_aligned(MSA,new_file):

	seq_num = len(MSA)
	length = len(MSA[-1])

	seqs =[]
	for x in MSA:
		if ('X' not in str(x.seq).upper()):
			new_seq = str(x.seq).upper()
			seqs.append(SeqRecord(Seq(new_seq,x.seq.alphabet), x.id,'',''))
	seq_num2 = len(seqs)
	MSA2= MultipleSeqAlignment(seqs)
	AlignIO.write(MSA2, "./"+new_file+"_data/"+new_file+"_Xed.fa", "fasta")
	logging.info(str("of " + str(seq_num) + " sequences starting, " + str(seq_num2) + " sequences were kept"))
	return seqs

	
def redundant(MSA, ident_matrix, threshold, new_file):

	length = MSA.get_alignment_length()
	seq_num = len(MSA)
	threshold = float(threshold)

	new_seqs = []
	used = {}
	to_remove = {}
	pairs =[]

	for x in MSA:
		used[x] = True
		if not(x in to_remove):
			redundant = {}
			for y in MSA:
				if (not(y in used) and not(y in to_remove)):
					if (ident_matrix[y.id.split("|")[0]][x.id.split("|")[0]] > threshold):
						pairs.append(y.id+" "+x.id+" %ident = "+str(ident_matrix[y.id.split("|")[0]][x.id.split("|")[0]])+". greater than threshold.")
						redundant[x.id] = True
						redundant[y.id] = True
			if len(redundant) >0:
				del redundant[random.choice(redundant.keys())]
				for z in redundant:
					to_remove[z] = True
	if len(pairs)>0:
		logging.info("the redundant sequenes are:"+"\n".join(pairs))
	for x in MSA:
		if not(to_remove.has_key(x.id)):
			new_seqs.append(SeqRecord(Seq(str(x.seq),x.seq.alphabet), str(x.id),'',''))				

	seq_num2 = len(new_seqs)
	MSA2= MultipleSeqAlignment(new_seqs)
	AlignIO.write(MSA2, "./"+new_file+"_data/global/"+new_file+"_non_redundant.fa", "fasta")
	logging.info(str("starting with: "+str(seq_num)+" and ending with: "+str(seq_num2)+" using a identity threshold: "+str(threshold)))
	plt.close()
	v = venn2(subsets=(seq_num,0,seq_num2),set_labels = ('input seq', '', 'non-redundant seq'))
	v.get_label_by_id("01").set_text("")
	v.get_patch_by_id('10').set_color('blue')
	v.get_patch_by_id('11').set_alpha(0.8)
	v.get_patch_by_id('11').set_color('purple')
	plt.savefig("./"+new_file+"_data/global/figures/redundancy_venn.png")
	plt.close()
	return MSA2

