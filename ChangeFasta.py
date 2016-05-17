# -*- coding: utf-8 -*-
"""
Created on Thu May  14 2016
Maxicircle Annotator
Compare a DNA sequences against the  genes in Maxicircle
@author: Said
"""
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse
import os
import glob

os.chdir("/Users/Said/Github/Maxicircle/Classification/DNA/Phylo/")

for file_in in glob.glob("*.genes.fna"):
	seq_record=SeqIO.parse(file_in,"fasta")	
	i=1
	fasta_out=[]
	for seq in seq_record:
		species=seq.description.split("|")[1]
		if len(species.split(" "))>1 :
			species=species.split(" ")[0][0]+"_"+species.split(" ")[1]
		description=seq.name+"_"+species+"_"+str(i)
		i=i+1
		seq.id=description
		seq.description=""
		fasta_out.append(seq)
		file_out=file_in.split(".")[0]+".fna"
		output_handle = open(file_out, "w")
		SeqIO.write(fasta_out, output_handle, "fasta")
		output_handle.close()

os.chdir("/Users/Said/Github/Maxicircle/Classification/AA/Phylo/")

for file_in in glob.glob("*.proteins.faa"):
	seq_record=SeqIO.parse(file_in,"fasta")	
	i=1
	fasta_out=[]
	for seq in seq_record:
		species=seq.description.split("|")[1]
		if len(species.split(" "))>1 :
			species=species.split(" ")[0][0]+"_"+species.split(" ")[1]
		description=file_in.split(".")[0]+"_"+species+"_"+str(i)
		i=i+1
		seq.id=description
		seq.description=""
		fasta_out.append(seq)
		file_out=file_in.split(".")[0]+".faa"
		output_handle = open(file_out, "w")
		SeqIO.write(fasta_out, output_handle, "fasta")
		output_handle.close()
