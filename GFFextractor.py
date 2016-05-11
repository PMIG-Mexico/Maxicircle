# -*- coding: utf-8 -*-
"""
Created on Wed May 11 12:31:55 2016
GFFextractor
@author: Said
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-reference",help="Fasta file with the Reference Chromosome, this version is built specifically for Maxicircle feature extraction",type=str)
parser.add_argument("-gff",help="GFF for feature extraction, this version is built specifically for Maxicircle feature extraction",type=str)
parser.add_argument("-out",help='Output filename',type=str)

args = parser.parse_args()

fasta_file=str(args.reference) #'Leishmania_mexicana_maxicircle.fasta'
in_file = str(args.gff) #'Lmex2.gff'
out_file= str(args.out) #'gffextractortest.fasta'

def reverseDNA(seq):
	'''
	Get the complementary strand of DNA
	'''
	reverse=[]
	DNA_complementary={'A':'T','T':'A','C':'G','G':'C','a':'t','t':'a','c':'g','g':'c'}
	for i in seq:
		try:
			reverse.append(DNA_complementary[i])
		except KeyError:
			reverse.append('n')
	out=''.join(reverse)[::-1]
	return(out)

class FastaExtractor:
	def __init__(self, in_file,fasta_file):
		fasta= SeqIO.parse(open(fasta_file),'fasta')
		refseq=str(fasta.next().seq)
		sequence_out=[]
		names_out=[]
		with open(in_file,'r')as f:
		    for line in f:
		    	if not line.startswith("#"):
		    		feature=line.split('\t')
		    		name=feature[0]
		    		start=int(feature[3])
		    		end=int(feature[4])
		    		strand=feature[6]
		    		if strand=='+':
		    			sequence_out.append(refseq[start:end])
		    			names_out.append(name)
		    		elif strand=='-':
		    			sequence_out.append(reverseDNA(refseq[start:end]))
		    			names_out.append(name)
		    		else:
		    			print('Feature '+ name+' with strand not defined. Options + or - ')
		self.sequences=sequence_out
		self.names=names_out

def WriteMultiFasta(multifasta,out_file):
	sequence_out=multifasta.sequences
	names_out=multifasta.names
	output_handle = open(out_file, "w")
	for i in range(len(sequence_out)):
		output_handle.write('>'+names_out[i])
		output_handle.write("\n")
		output_handle.write(sequence_out[i])
		output_handle.write("\n")
	output_handle.close()
	print(str(i+1) +' sequences extracted')


multifasta=FastaExtractor(in_file, fasta_file)
WriteMultiFasta(multifasta,out_file)
