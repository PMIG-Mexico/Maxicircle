# -*- coding: utf-8 -*-
"""
Created on Thu May  25 2016

@author: Said
"""
from Bio import SearchIO
import argparse
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbiblastxCommandline
import glob
import csv


class BlastMaxicircle:
	def __init__(self,Maxicircle,out_file):
		file_out=open(out_file,'w')
		writer = csv.writer(file_out,delimiter="\t")
		writer.writerow(["##gff-version","3"])
		rows=[]
		
		for protein in glob.glob("/Users/Said/Github/Maxicircle/DB/AA/*.faa"):
			output_file=protein.split("/")[-1]+".xml"
			blastx_cline = NcbiblastxCommandline(query=Maxicircle , db=protein, 
		                                      outfmt=5, out=output_file)
			blastx_cline()
			blast_qresult = SearchIO.read(output_file, 'blast-xml')
			if len(blast_qresult)>0:
				best=blast_qresult[0][0]
				query_range=[x for x in best.query_range]
				if best.query_strand>0:
					query_strand="+"
				else:
					query_strand="-"
				chromosome=best.query_id
				rows.append([chromosome,".","exon",query_range[0],query_range[1],".",query_strand,".","ID="+protein.split("/")[-1].split(".faa")[0]])
		
		print(str(len(rows))+" exons found")
		rows=iter(rows)
		writer.writerows(rows)
		
def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("-fasta",help="Maxicircle Fasta",type=str,required=True)
	parser.add_argument("-out_file",help='Output filename',type=str,default='Maxicircle.gff')
	args = parser.parse_args()
	Maxicircle_file=args.fasta
	out_file=args.out_file
	BlastMaxicircle(Maxicircle_file, out_file)
	
if __name__ == "__main__":
	main()

