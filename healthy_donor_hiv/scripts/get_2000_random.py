#!/usr/bin/env python
from Bio import SeqIO
import random
import sys
handle = open(sys.argv[1])
records = SeqIO.parse(handle, "fasta")
dictionary = {}
for seq in records:
   dictionary[seq.id] = str(seq.seq)
random_dict = random.sample(dictionary.items(),2000)

with open(sys.argv[2],'w') as f:
    for seq in random_dict:
        f.write(">{0}\n{1}\n".format(seq[0],seq[1]))
