from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
from Bio import SeqIO
import urllib,urllib2
import csv
import logging
import time
import ETA
import gzip

def uniprot(txt):
	return txt.split('_')[0]

def Pfam(input_file, run_name):
	seqs =[]
	of = len(input_file)
	on = 0
	time0 = time.time()
	for x in input_file:
		on = on +1
		print "Getting species: "+ETA.ETA(on,of,time0,"N")
		url = 'http://www.uniprot.org/uniprot/'+uniprot(x.id)+'.xml'	
		time.sleep(2)
		try:
			response = urllib2.urlopen(url,timeout=900)
		except urllib2.URLError, e:
			logging.info("error retrieving: "+uniprot(x.id))
		except:
			logging.info("trouble downloading: "+uniprot(x.id))
		else:
			file = open('./'+run_name+'_data/xml/'+x.id+'.xml', "w")
			file.write(response.read())
			file.close()
			file = open('./'+run_name+'_data/xml/'+x.id+'.xml', 'rb')
			try:
				new_seq = SeqIO.read(file, "uniprot-xml")
			except:
				logging.info("cannot load xml file for: "+uniprot(x.id))
			else:
				try:
					new_spec = new_seq.annotations['organism'] 
				except:
					logging.info("no species info for: "+uniprot(x.id))
				else:
					if not('sequence_fragment' in new_seq.annotations):		
						b = '_'.join(new_spec.split()[0:2])
						logging.info("got record for:"+str(uniprot(x.id)+"|"+b))
						seqs.append(SeqRecord(Seq(str(x.seq),x.seq.alphabet), str(uniprot(x.id)+"|"+b),'',''))
					else:
						logging.info("fragment, skipping: "+uniprot(x.id))
	MSA2= MultipleSeqAlignment(seqs)
	AlignIO.write(MSA2, "./"+run_name+"_data/"+run_name+"_raw_seq.fa", "fasta")
	return seqs	

def Pfam_unaln(input_file, run_name):
	seqs =[]
	of = len(input_file)
	on = 0
	time0 = time.time()
	for x in input_file:
		on = on +1
		print "Getting species: "+ETA.ETA(on,of,time0,"N")
		url = 'http://www.uniprot.org/uniprot/'+uniprot(x.id)+'.xml'	
		time.sleep(2)
		try:
			response = urllib2.urlopen(url,timeout=900)
		except urllib2.URLError, e:
			logging.info("error retrieving: "+uniprot(x.id))
		except:
			logging.info("trouble downloading: "+uniprot(x.id))
		else:
			file = open('./'+run_name+'_data/xml/'+x.id+'.xml', "w")
			file.write(response.read())
			file.close()
			file = open('./'+run_name+'_data/xml/'+x.id+'.xml', 'rb')
			try:
				new_seq = SeqIO.read(file, "uniprot-xml")
			except:
				logging.info("cannot load xml file for: "+uniprot(x.id))
			else:
				try:
					new_spec = new_seq.annotations['organism'] 
				except:
					logging.info("no species info for: "+uniprot(x.id))
				else:
					if not('sequence_fragment' in new_seq.annotations):		
						b = '_'.join(new_spec.split()[0:2])
						logging.info("got record for:"+str(uniprot(x.id)+"|"+b))
						seqs.append(SeqRecord(Seq(str(x.seq),x.seq.alphabet), str(uniprot(x.id)+"|"+b),'',''))
					else:
						logging.info("fragment, skipping: "+uniprot(x.id))
	SeqIO.write(seqs, "./"+run_name+"_data/"+run_name+"_raw_seq.fa", "fasta")
	return seqs	
