#!/usr/bin/env python
from optparse import OptionParser
import os
import sys
import pipes


usage = "set_pdb_chain_id.py original_chain new_chain input_file.pdb output_file.pdb"

if len(sys.argv) != 5:
    print("specify original chain ID, desired new chain ID, input file name, and output file name")
    print(usage)
old_chain = sys.argv[1]
new_chain = sys.argv[2]
pdb_file = open(sys.argv[3],'r')
out_file = open(sys.argv[4],'w')

for line in pdb_file:
    if line[0:4] != "ATOM":
        continue
    line = list(line)
    if(line[21] == old_chain):
        line[21] = new_chain
    line = "".join(line)
    out_file.write(line)
pdb_file.close()
out_file.close()
os.system("/sb/meiler/scripts/capture_command.sh " + ' '.join([pipes.quote(x) for x in sys.argv]))


