#!/usr/bin/env python
from Bio import SeqIO as so

import sys

l = [ i for i in so.parse(sys.argv[1],'fasta')]



print l[0]
for i in l:
	if i.name != 'PG9' and i.name != 'PG16':
		#i.id = ':'.join(i.id.split(':')[5:])
		i.id = i.name
		print i.name, i.id
	else:
	    print "print found", i.id
            continue

start = "use_input_sc\nex 1 ex 2\nNATAA\nstart\n"
for i in l:
	handle=open(i.id+'.resfile','w').write("")
	handle=open(i.id+'.resfile','a')
	handle.writelines(start)
	count = 1
	for seq in i.seq:
		if ((count == 110 or count == 111)and seq == 'Y'):
			handle.write(str(count)+"\t"+"H\t"+"NATAA\n")
			count += 1
		else:
		   handle.write(str(count)+"\t"+"H\t"+"PIKAA\t"+seq+"\n")
		   count +=1

