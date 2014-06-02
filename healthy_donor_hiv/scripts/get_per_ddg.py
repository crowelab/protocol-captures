#!/usr/bin/env python
from Bio.PDB import PDBParser as pdbparser
from amino_acids import longer_names
import warnings
import argparse
import os
from Bio.PDB import PDBExceptions
from multiprocessing import Pool
warnings.simplefilter('ignore', PDBExceptions.PDBConstructionWarning)

# usage
usage = "%prog [options] -n native.pdb *pdbs"
parser = argparse.ArgumentParser(prog="get_per_ddg.py", description="The Per Residue DDG mover outputs a poorly parseble result at the end of pdbs.\
                                                                    Use this script to take in all the pdb files and output a nice tab seperated file\
                                                                    of each residue number, chain, and identity along with it's per residue ddg. It \
                                                                    will also grab the total ddg for you and store it another file. Each pdb file you pass\
                                                                    is given as an id that you can look up in the other file. This is so the output doesn't get cluttered\
                                                                    trust me, you will like it\n\
                                                                    THIS SCRIPT ASSUMES YOU HAVE RUN BOTH THE DDG FILTER AND THE PER_DDG_MOVER FROM ROSETTA SCRIPTS!!!!")

parser.add_argument("pdbs", metavar="*pdb", nargs=argparse.REMAINDER, help="The pdb files that you want superimposed")
parser.add_argument("--out", "-o", dest="table", help="the name of the prefix of the  output files", default="ddgs")
parser.add_argument("--multi", "-m", dest="multi", help="run using multiple processors", action='store_true', default=False)
args = parser.parse_args()

if len(args.pdbs) < 1:
    parser.error("specify at least 1 protein to extract ddgs")

files = args.pdbs


def get_ddg_table(file):
    dict_of_ddgs = {}
    for line in open(file).readlines():

        try:
            splitter = line.split()[0].split("_")
            if splitter[0] == "residue":
                pose_number = splitter[2]
                energy = line.split()[1]
                dict_of_ddgs[int(pose_number)] = energy
        except KeyError:
            continue

        try:
            if line.split()[0] == "ddg":
                total_ddg = line.split()[1]
        except KeyError:
            continue
    return dict_of_ddgs, total_ddg

parser = pdbparser()


def get_ddg(pdb):
    returner_total_ddg = {}
    returner_per_ddg = {}
    pdb_file = parser.get_structure(pdb, pdb)
    residues = [x for x in pdb_file.get_residues()]
    ddg_table, total_ddg = get_ddg_table(pdb)
    returner_total_ddg[os.path.basename(pdb)] = total_ddg
    #super_out.write("{0}\t{1}\n".format(os.path.basename(pdb), total_ddg))
    for index, residue in enumerate(residues, start=1):
        chain = residue.get_full_id()[2]
        res_number = residue.get_full_id()[-1][1]
        pose_number = index
        try:
            residue_id = longer_names[residue.resname]
        except KeyError:
            residue_id = residue.resname
        #out_file.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(
        #    os.path.basename(pdb), chain, res_number, pose_number, residue_id, ddg_table[pose_number]))
        returner_per_ddg[(os.path.basename(pdb),chain,res_number,pose_number,residue_id)] = ddg_table[pose_number]
    return returner_total_ddg, returner_per_ddg

def unwrapper(pdb_super_out):
    return get_ddg(*pdb_super_out)


def write_out_files(total_ddg,per_ddg,out_file,super_out):
    super_out.write("{0}\t{1}\n".format(total_ddg.keys()[0],total_ddg.values()[0]))
    for residue in sorted(per_ddg):
        model = residue[0]
        chain = residue[1]
        res_num = residue[2]
        pose_num = residue[3]
        identity = residue[4]
        ddg = per_ddg[residue]
        out_file.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(model,chain,res_num,pose_num,identity,ddg))


def main():
    if args.multi:
        try:
            p = Pool()
            ddg = p.map(get_ddg, files)
        except KeyboardInterrupt:
            p.terminate()
            quit()
        with open("total"+args.table+".tsv", 'w') as super_out:
            super_out.write("{0}\t{1}\n".format("Model", "total_ddg"))
            with open("per_residue_"+args.table+".tsv",'w') as out_file:
                out_file.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format("Model", "chain", "residue_num", "pose_number", "residue_id", "per_ddg"))
                for index in ddg:
                    total_ddg = index[0]
                    per_ddg = index[1]
                    write_out_files(total_ddg,per_ddg,out_file,super_out)
    else:
        with open("total"+args.table+".tsv", 'w') as super_out:
            super_out.write("{0}\t{1}\n".format("Model", "total_ddg"))
            with open("per_residue_"+args.table+".tsv",'w') as out_file:
                out_file.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format("Model", "chain", "residue_num", "pose_number", "residue_id", "per_ddg"))
                for file in files:
                    total_ddg, per_ddg = get_ddg(file)
                    super_out.write("{0}\t{1}\n".format(total_ddg.keys()[0],total_ddg.values()[0]))
                    for residue in sorted(per_ddg):
                        model = residue[0]
                        chain = residue[1]
                        res_num = residue[2]
                        pose_num = residue[3]
                        identity = residue[4]
                        ddg = per_ddg[residue]
                        out_file.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(model,chain,res_num,pose_num,identity,ddg))

if __name__ == '__main__':
    main()
    os.system("/sb/meiler/scripts/capture_command.sh " + ' '.join([pipes.quote(x) for x in sys.argv]))
