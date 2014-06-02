#!/usr/bin/env python2.6

import sys
from Bio import SeqIO as IO

if (len(sys.argv)<2):
	print "Usage: script <catenated fasta_file> <replacement start> <replacement end> <chain> <left_flank_lenght_length>\ncatenated fasta_file=the fastafiles you want to turn into individual residue files\nreplacement start=the pdb res id where you want to start threading\nchain=the chain you want to thread over in your PDB\n"
	exit()
Fastahandle = sys.argv[1]


for seq_record in IO.parse(open(Fastahandle), "fasta"):
	#numbers = [i for i in range(int(sys.argv[2]),int(sys.argv[3]))]
	#print numbers
#	count = int(sys.argv[2])
	open(seq_record.description + ".resfile", 'w').write("#"+ str(seq_record.description) + "\nUSE_INPUT_SC\nNATAA\nEX 1 EX 2\nstart\n")
	for num,aa in enumerate(seq_record.seq, start=96):
		#if num >= int(sys.argv[5]) and count <= int(sys.argv[3]):
			#if (count == 111 or count == 110) and aa == "Y":
			#	design = " NATAA "
			#else:
			design = " PIKAA " + str(aa)
			open(seq_record.description + ".resfile", 'a').write(str(num) + "\t" + "H" + design +"\n")
		#	count+=1
