#!/usr/bin/env python2.7
import sys, os, pipes
import rosettaScore
from Bio.PDB import *
from optparse import OptionParser
import warnings
def main():
    usage = "%prog input.pdb output.pdb"
    parser= OptionParser(usage)
    parser.add_option("-n",dest="start",help="residue number to start with, default is 1",default=1)
    parser.add_option("--preserve",dest="preserve",help="preserve insertion code and heteroflags",default=False, action="store_true")
    parser.add_option("--norestart",dest="norestart",help="don't start renumbering at each chain, default=False",default=False, action="store_true")
    parser.add_option("--keep-table",dest="table",help="Preserve the rosetta score table at the bottom of the pdb",action="store_true",default=False)
    (options, args) = parser.parse_args()

    if len(args) < 1:
        parser.error("Arguments required. Use -h or --help")

    warnings.simplefilter('ignore',PDBExceptions.PDBConstructionWarning)
    #cores = rosettaScore.ScoreTable(args[0])

    PDBparse = PDBParser(PERMISSIVE=1)
    struct = PDBparse.get_structure(args[0][0:3],args[0])
    residue_id = int(options.start)
    chain_id = ""
    for residue in struct.get_residues():
        chain = residue.get_parent()
        if(chain_id != chain.get_id() and not options.norestart):
            chain_id = chain.get_id()
            residue_id=int(options.start)
        #print chain.get_id()
        if(options.preserve):
            hetero = residue.id[0]
            insert = residue.id[2]
            residue.id=(hetero,residue_id,insert)
        else:
            residue.id=(' ',residue_id,' ')
        residue_id +=1


    io=PDBIO()
    io.set_structure(struct)
    io.save(args[1])

    if(options.table):
        raw_table = rosettaScore.get_table(args[0])
        outfile = open(args[1],'a')
        outfile.writelines(raw_table)
        outfile.close()


if __name__ == "__main__":
    main()
    # add the command line arguments to the database
    os.system("/sb/meiler/scripts/capture_command.sh " + ' '.join([pipes.quote(x) for x in sys.argv]))

