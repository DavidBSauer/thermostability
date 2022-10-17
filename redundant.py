import random
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
import logging
import ETA

def remover(MSA, ident_matrix, threshold, new_file):

	length = MSA.get_alignment_length()
	seq_num = len(MSA)
	threshold = float(threshold)

	new_seqs = []
	used = {}
	to_remove = {}
	pairs =[]
	count =0
	for x in MSA:
		print "Removing redundant sequences: "+ETA.ETA(count,seq_num,time0,"1/2N^2")
		used[x] = True
		if not(x in to_remove):
			redundant = {}
			for y in MSA:
				if (not(y in used) and not(y in to_remove)):
					if (ident_matrix[y.id][x.id] > threshold):
						pairs.append(y.id+" and "+x.id+" %ident = "+str(ident_matrix[y.id][x.id])+". greater than threshold.")
						redundant[x.id] = True
						redundant[y.id] = True
			if len(redundant) >0:
				del redundant[random.choice(redundant.keys())]
				for z in redundant:
					to_remove[z] = True
	logging.info("\n".join(pairs))
	for x in MSA:
		if not(to_remove.has_key(x.id)):
			new_seqs.append(SeqRecord(Seq(str(x.seq),x.seq.alphabet), str(x.id),'',''))				

	seq_num2 = len(new_seqs)
	MSA2= MultipleSeqAlignment(new_seqs)
	AlignIO.write(MSA2, "./"+new_file+"_data/"+new_file+"_non_redundant.fa", "fasta")
	logging.info("starting with: "+str(seq_num)+" and ending with: "+str(seq_num2)+" using a identity threshold: "+str(threshold))
	return MSA2
