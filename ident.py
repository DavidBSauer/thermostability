# ident_calculator

from Bio.Align import MultipleSeqAlignment
import logging
import ETA
import time



def calc(MSA,out_file):
	logging.info("calculating ident matrix")
	tried = {}
	matrix = {}	
	count = 0
	num = len(MSA)
	time0 = time.time()
	for x in MSA:
		count = count+1
		print "Calculating identity matrix: "+ETA.ETA(count,num,time0,"1/2N^2")		
		line = {}
		for y in MSA:
			if not(y.id in tried):
				#calc percent ident
				length = 0
				ident = 0
				for z in range(0,len(y.seq)):
					if ((x.seq[z] != '-') or (y.seq[z] != '-')):
						length = length +1
						if (x.seq[z] == y.seq[z]):
							ident = ident+1				
				line[y.id.split("|")[0]] = 100*float(ident)/float(length)
			else:
				line[y.id.split("|")[0]] = matrix[y.id.split("|")[0]][x.id.split("|")[0]]
		tried[x.id] = True
		matrix[x.id.split("|")[0]] = line
		
	f = open("./"+out_file+"_data/"+out_file+"_ident_matrix.txt","w+")	
	lokeys = matrix.keys()
	for x in lokeys:
		f.write(str(x)+"\t")
		for y in range(0,len(lokeys)):
			f.write(str(matrix[x][lokeys[y]])+"\t")
		f.write("\n")
	f.close()
	logging.info("have ident matrix for "+str(len(matrix.keys()))+" sequences")

	return matrix

''''
	A	B	C
A	0	1	2
B	1	0	3
C	2	3	0

pass 1:
A:A - 0
A:B - 1
A:C - 2
A= tried
pass 2:
B:A - tried -> 1
B:B - 0
B:C - 3
B = tried
pass 3
C:A - tried -> 2
C:B - tried -> 3
C:C -0

file:
A	0	1	2
B	1	0	3
C	2	3	0

'''		

def calc_already_done(file):
	logging.info("using precomputed ident matrix")
	f = open(file,"r")
	matrix = {}
	lokeys = []
	for work_line in f:
		lokeys.append(work_line.split()[0])
	f.close()
	f = open(file,"r")
	for work_line in f:
		split_line = work_line.split()
		first = split_line.pop(0)
		line={}
		for x in range(0,len(split_line)):
			line[lokeys[x]]=float(split_line[x])
		matrix[first] = line
	f.close()
	logging.info("have ident matrix for "+str(len(matrix.keys()))+" sequences")
	return matrix
