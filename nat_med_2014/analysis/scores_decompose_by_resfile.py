#!/usr/bin/env python2.7
import sys
import os
import warnings
from multiprocessing import Pool
import argparse
import amino_acids
import rosettaScore_beta as rsb
from Bio.PDB import PDBExceptions
from rosettautil import resfile

warnings.simplefilter('ignore', PDBExceptions.PDBConstructionWarning)

if sys.version_info < (2, 7):
    raise Exception("You must use python2.7 to run this")


# usage
usage = "%prog [options] -n native.pdb *pdbs"
parser = argparse.ArgumentParser(prog="score_decompose_by_resfile.py", description="get scores for every residue defined in the resfile")

parser.add_argument("pdbs", metavar="*pdb", nargs=argparse.REMAINDER,help="The pdb files that you want superimposed")
parser.add_argument(
    "--multi", "-m", dest="multi", help="run using multiple processors", action='store_true', default=False)
parser.add_argument("--out", "-o", dest="table", help="the output table", default="output_scores.tsv")
parser.add_argument("--res", "-r", dest="residues",required=True,
                    help="the res file to align by too, make sure this resfile corresponds to the decoy. If not specified, it will align the whole molecule", default="")
args = parser.parse_args()



# get the clustal correlations to be used in the rmsd calculations
#if len(args.pdbs) == 1:
 #   args.pdbs = [args.pdbs]

def main(args):
    if args.multi:
        try:
            p = Pool()
            tag_list = p.map(find_scores, args.pdbs)
        except KeyboardInterrupt:
            p.terminate()
            quit()
    else:
        tag_list = []
        for i in args.pdbs:
            # stats is a dictionary
            stats = find_scores(i)
            tag_list.append(stats)
    make_table(tag_list, args.table)


def make_table(tag_list, out_name):
    additional_headers = tag_list[0].values()[0].keys()
    header = ["MODEL", "CHAIN", "RES_NUM", "RESIDUE"] + additional_headers + ["\n"]
    with open(out_name,'w') as f:
        f.write("\t".join(header))
        for entity in tag_list:
            for model_entity in entity:
                model = os.path.basename(model_entity[0])
                chain = model_entity[1]
                res_num = model_entity[2]
                residue = model_entity[3]
                f.write(model + "\t")
                f.write("{0}\t".format(chain))
                f.write("{0}\t".format(res_num))
                f.write("{0}\t".format(residue))
                for scores in entity[model_entity]:
                    f.write("{0:.2f}\t".format(entity[model_entity][scores]))
                f.write("\n")


def find_scores(model):
    residue_scores = {}
    score_table = rsb.ScoreTable(model)
    residues = resfile.Residue_File(args.residues).get_designed_entities()
    for residue in residues:
        chain = residue[0]
        pdbres = int(residue[1])
        score_table_per_residue = score_table.get_all_score_terms(chain=chain, pdbres=pdbres)
        scores = {}
        for term in score_table_per_residue.iterkeys():
            look_up = score_table.get_score(term=term,chain=residue[0],pdbres=int(residue[1]))
            score = look_up[0]
            identity = amino_acids.longer_names[look_up[1]]
            scores[term] = score
        key = (model,chain,pdbres,identity)
        residue_scores[key] = scores
    return residue_scores


if __name__ == "__main__":
    main(args)
