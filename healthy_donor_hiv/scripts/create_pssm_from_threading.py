#!/usr/bin/env python2.7
import rosettaScore_beta as rs
import sys
from amino_acids import longer_names
import argparse
import textwrap

def args():
		parser = argparse.ArgumentParser(prog="create_pssm_from_threading.py",formatter_class=argparse.RawTextHelpFormatter,description=textwrap.dedent('''\
			Create PSSM
		------------------------------------------------------------\n
		If you have a bunch of PDBs you want to make a PSSM out of their energies, simply use this script. 
		It will generate a Position Structure Specific Scoring Matrix (P3SM). 
		It can then read in a P3SM along with raw fasta files to predict what their score would be given the training of the P3SM.
		Author - Jordan Willis'''))
	
		neces = parser.add_argument_group(title='Necessary',description="These have to be included")
		neces.add_argument("targets",metavar="\*.pdbs or \*.fasta", nargs=argparse.REMAINDER,help="If you are in generation mode, you need to include the PDBs that should be scored. \
		If you are in scoring mode, you must include a fasta file including all the fasta file syou want scored")
		neces.add_argument("-s",dest="score_mode",action="store_true",help="score a bunch of fasta sequences given a P3SM")
		neces.add_argument("-g",dest="generate_mode",action="store_true",help="generate a p3sm given a bunch of pdbs")
		neces.add_argument("-p",dest="p3sm",help="the p3sm you want to use if you are in score mode to score all the fastas")
		output_ops = parser.add_argument_group(title="optionL", description="optional but really recommended")
		output_ops.add_argument("-r","--resfile",dest="resfile",help="The resfile you want to score by or just generate your p3sm. makes it so much easier.")
		output_ops.add_argument("-n","--name",dest="name",default="output.p3sm",help="the name of the p3sm you want output in generation mode")
		return parser.parse_args()


class p3sm():
	def _parse_resis(self,res):
		res_list = []
		try:
			with open(res) as f:
				breaker = False
				for z,line in enumerate(f.readlines()):
					if line.split()[0] =="start":
						breaker = True
						continue
					if breaker:
						entity = z
						res_number = int(line.split()[0])
						chain = line.split()[1]
						res_list.append((chain,res_number))
				return res_list
		except IOError:
			print "res file didn't parse"
			return ""
		
	
	def __init__(self,resids="",term="total"):
			self.resis = self._parse_resis(resids)
			self.term = term
			self.create_matrix = {}
			self.count_matrix = {} #for making an average
	
	def generate_pssm(self,files):
			self.create_pssm(files)
			self.average_pssm()
					

	
	def create_pssm(self,files):
		counter = 0
		for l,i in enumerate(files):
			scores = rs.ScoreTable(i)
			for j in self.resis:
				current_id = scores.get_score(chain=j[0],term=self.term,pdbres=j[1])
				position = j[1]
				#print position, current_id[1]
				name = longer_names[current_id[1]]
				position_score = float(current_id[0])
				try:
					self.create_matrix[position][name] += position_score
					self.count_matrix[position][name] += 1
				except KeyError:
					self.create_matrix[position] = {'A':0,'R':0,'N':0,'D':0,'C':0,'E':0,'Q':0,'G':0,'H':0,'I':0,'L':0,'K':0,'M':0,'F':0,'P':0,'S':0,'T':0,'W':0,'Y':0,'V':0}
					self.count_matrix[position] = {'A':0,'R':0,'N':0,'D':0,'C':0,'E':0,'Q':0,'G':0,'H':0,'I':0,'L':0,'K':0,'M':0,'F':0,'P':0,'S':0,'T':0,'W':0,'Y':0,'V':0}
					self.create_matrix[position][name] += position_score
					self.count_matrix[position][name] += 1
			if l % 100 == 0:
				counter += 100
				print "{0} files done in generating pssm".format(counter)
					
	def average_pssm(self):
		for i in self.resis:
			for j in self.create_matrix[i[1]]:
				try:
					self.create_matrix[i[1]][j] = self.create_matrix[i[1]][j]/self.count_matrix[i[1]][j]
				except ZeroDivisionError:
					self.create_matrix[i[1]][j] = 0

	def output_pssm(self,filename="outfile_pssm.psssm"):
		positions = self.create_matrix.keys()
		aa = self.create_matrix[positions[0]].keys()
		aa = ["{0:^}".format(x) for x in aa]
		with open(filename,'w') as f:
			f.write("\t  "+"   ".join(aa)+"\n")
			for i in positions:
				f.write("{0}\t".format(i))
				row = []
				for j in aa:
					row.append("{0:.2f}".format(self.create_matrix[i][j]))
				f.write(" ".join(row)+"\n")
	
	def read_in_pssm(self,input):
		file = open(input).readlines()
		aa = [x for x in file[0].split() if x !=""]
		for i in file[1:]:
			line = i.split()
			position = line[0]
			for a,j in zip(aa,line[1:]):
				try:
					self.create_matrix[position][a] = j
				except:
					self.create_matrix[position] = {'A':0,'R':0,'N':0,'D':0,'C':0,'E':0,'Q':0,'G':0,'H':0,'I':0,'L':0,'K':0,'M':0,'F':0,'P':0,'S':0,'T':0,'W':0,'Y':0,'V':0}
					self.create_matrix[position][a] = j
					
	
	def score_fasta_files(self,input,output):
		from Bio import SeqIO as so
		file = so.parse(input,'fasta')
		with open(output,'w') as f:
			for fasta in file:
				id = fasta.id
				sequence = fasta.seq
				score = 0.0
				for position,letter in enumerate(sequence,start=self.resis[0][1]):	
					score += float(self.create_matrix[str(position)][letter])
				f.write("{0}\t{1}\t{2}\n".format(id,sequence,score))
	
	def debug(self):
		for i in self.create_matrix:#sorted([int(x) for x in self.create_matrix.iterkeys()]):
				print i, self.create_matrix[i]	



			
if __name__ == '__main__':
	args = args()
	c = p3sm(resids=args.resfile)
	print "resfile_parsed"
	if args.generate_mode:
		c.generate_pssm(args.targets)
		c.output_pssm(filename=args.name)
	elif args.score_mode:
		c.read_in_pssm(input=args.p3sm)
		c.score_fasta_files(input=args.targets,output="scored_fasta.output")
	else:
		print "you need a score mode or generate mode use -h"
