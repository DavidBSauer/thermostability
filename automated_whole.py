## automated method for doing entire global analysis
#save files each step but doing whole thing automatically except when there are exceptions

import sys
import csv

from sys import argv
import os
import getter
import sorter
import matrix_stat as matrix
from Bio import AlignIO
from Bio import SeqIO
import ident
import remover
import pairwise
import time
import logging
import datetime

def main():
	script, MSA_file, ref_file, run_name, input_type = argv
	t0 = time.time()
	folder = "./"+run_name+"_data"
	if not os.path.exists(folder):
 	   os.makedirs(folder)

	logging.basicConfig(filename=str(folder+"/"+run_name+'.log'), level=logging.INFO)
 

	#print basic parameters
	logging.info("starting at: "+time.strftime("%c"))
	logging.info("input file: "+MSA_file)
	logging.info("ref file: "+ref_file)
	logging.info("run name: "+run_name)
	logging.info("input type: "+input_type)
	if not os.path.exists(folder+"/xml/"):
		os.makedirs(folder+"/xml/")

	
	if input_type != "SAO":
		print "Give the analysis type (Global: G, Pairwise: P, or Both: B):"
		mode = raw_input().upper()
		logging.info("mode: "+mode)
		if (mode == "G"):
			print "Give the redundancy treshold for the global calculation:"
			global_ident = float(raw_input())
			logging.info("global ident threshold: "+str(global_ident))
			if not os.path.exists(folder+"/global/"):
				os.makedirs(folder+"/global/")
			if not os.path.exists(folder+"/global/figures/"):
				os.makedirs(folder+"/global/figures/")
			if not os.path.exists(folder+"/global/figures/histograms/"):
				os.makedirs(folder+"/global/figures/histograms/")
		elif (mode == "P"):
			print "Give the temperature treshold for the pairwise calculation:"
			pair_temp = float(raw_input())
			logging.info("pairwise temp threshold: "+str(pair_temp))
			print "Give the identity treshold for the pairwise calculation:"
			pair_ident = float(raw_input())
			logging.info("pairwise ident threshold: "+str(pair_ident))
			print "Give the gap frequency treshold for the pairwise calculation:"
			gap = float(raw_input())
			logging.info("pairwise gap threshold: "+str(gap))
			if not os.path.exists(folder+"/pairwise/"):
				os.makedirs(folder+"/pairwise/")
			if not os.path.exists(folder+"/pairwise/figures/"):
				os.makedirs(folder+"/pairwise/figures/")
			if not os.path.exists(folder+"/pairwise/figures/histograms/"):
				os.makedirs(folder+"/pairwise/figures/histograms/")
		elif (mode == "B"):
			print "Give the redundancy treshold for the global calculation:"
			global_ident = float(raw_input())
			logging.info("global ident threshold: "+str(global_ident))
			print "Give the temperature treshold for the pairwise calculation:"
			pair_temp = float(raw_input())
			logging.info("pairwise temp threshold: "+str(pair_temp))
			print "Give the identity treshold for the pairwise calculation:"
			pair_ident = float(raw_input())
			logging.info("pairwise ident threshold: "+str(pair_ident))	
			print "Give the gap frequency treshold for the pairwise calculation:"
			gap = float(raw_input())
			logging.info("pairwise gap threshold: "+str(gap))
			if not os.path.exists(folder+"/global/"):
				os.makedirs(folder+"/global/")
			if not os.path.exists(folder+"/global/figures/"):
				os.makedirs(folder+"/global/figures/")
			if not os.path.exists(folder+"/global/figures/histograms/"):
				os.makedirs(folder+"/global/figures/histograms/")
				if not os.path.exists(folder+"/pairwise/"):
					os.makedirs(folder+"/pairwise/")
				if not os.path.exists(folder+"/pairwise/figures/"):
					os.makedirs(folder+"/pairwise/figures/")
				if not os.path.exists(folder+"/pairwise/figures/histograms/"):
					os.makedirs(folder+"/pairwise/figures/histograms/")
		else:
			print "give a valid analysis format"
			logging.info("did not give a valid analysis format")
			sys.exit()

	logging.info("trying to open species-temp file")
	infile = open(ref_file,'r')
	reader = csv.reader(infile,delimiter='\t')
	species_temp = dict((str(rows[0]),float(rows[1])) for rows in reader)
	infile.close()    
	logging.info("found "+str(len(species_temp.keys()))+" species-OGT pairs")

	
	if (input_type == "SAO"):
		logging.info("opening sequence file")
		MSA = AlignIO.read(MSA_file,"fasta")
		logging.info("sequence file had "+str(len(MSA))+" sequences")
		#use a provided alignment file, just get and keep those with assignment. in Fasta format
		logging.info("done\nusing method where just takes file and retrieves species. not full calculation\ngetting sequences")
		seqs = getter.Pfam_unaln(MSA,run_name)
		logging.info("done\nremoving Xs/gaps & convert to upper case...")
		seqs = remover.cleanup(seqs,run_name)
		logging.info("done\nkeeping only those with temp data...")
		sorter.assign(seqs,species_temp,run_name)
		logging.info("done\nexiting.")
		sys.exit()		
	elif (input_type == "ALN"):
		logging.info("opening MSA file")
		MSA = AlignIO.read(MSA_file,"fasta")
		logging.info("MSA file had "+str(len(MSA))+" sequences")
		logging.info("done\nusing method where an alignment is already provided.\ngetting species")
		#use a provided alignment file, in Fasta format
		seqs = getter.Pfam(MSA, run_name)
		logging.info("done\nemoving Xs & convert to upper case...")
		seqs = remover.cleanup_already_aligned(seqs,run_name)
		logging.info("done\nassigning temps and sorting...")
		seqs = sorter.sort(seqs,species_temp,run_name)
		logging.info("done\ncalc ident matrix...")
		idents = ident.calc(seqs,run_name)
	elif (input_type == "STA"):
		logging.info("opening MSA file")
		MSA = AlignIO.read(MSA_file,"fasta")
		logging.info("MSA file had "+str(len(MSA))+" sequences")
		logging.info("done\nusing method where an asigned alignment and ident matrix are already provided.\ncalc ident matrix...")
		idents = ident.calc_already_done(ref_file)
		seqs = MSA
	else:
		print "give a valid data input format"
		logging.info("did not give valid input format")
		sys.exit()
	
	if (mode == "G"):
		#running in global mode (separating into thermophiles/mesophiles and comparing AA frequency)
		logging.info("done\nremove redundant sequences")
		seqs = remover.redundant(seqs,idents,global_ident,run_name)
		logging.info("done\nmatrix comparison...")
		matrix.divide(seqs,run_name)		
		logging.info("done\nexiting.")
	
	elif (mode == "P"):
		#running in pairwise mode (comparing sequences pairwise and calc frequency of positional differences)
		logging.info("done\npairwise comparison...")
		pairwise.analyse(seqs,idents,pair_temp,pair_ident,gap,run_name)
		logging.info("done\nexiting.")

	elif (mode == "B"):
		#doing both pairwise and global comparisons
		logging.info("done\npairwise comparison...")
		pairwise.analyse(seqs,idents,pair_temp,pair_ident,gap,run_name)
		
		#do global analysis second
		logging.info("done\nremove redundant sequences")
		seqs = remover.redundant(seqs,idents,global_ident,run_name)
		logging.info("done\nmatrix comparison...")
		matrix.divide(seqs,run_name)
		logging.info("done\nexiting.")

	logging.info("finished normally at: "+time.strftime("%c"))
	logging.info("total time: "+str(datetime.timedelta(seconds=time.time()-t0)))

if __name__ == '__main__':
    main()
