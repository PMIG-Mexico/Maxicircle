# -*- coding: utf-8 -*-
"""
Created on Thu May  5 16:24:25 2016

@author: Said
"""
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import ProteinAlphabet
import argparse

class ProteinSeqFromGenBank:
	def __init__(self,GeneBank_file):
		fasta_out=[]
		for seq_record in SeqIO.parse(GeneBank_file,"gb"):
			name="".join(seq_record.name)
			organism="".join(seq_record.annotations['organism'])
			for feature in seq_record.features:
				if 'translation' in feature.qualifiers:
					proteinSequence=feature.qualifiers['translation']
					proteinID=feature.qualifiers['protein_id']
					if 'product' in feature.qualifiers:
						product=feature.qualifiers['product']
					else:
						product=['']*len(proteinSequence)
					for i in range(len(proteinSequence)):
						description='|'+organism+'|'+product[i]
						out=SeqRecord(seq=Seq(proteinSequence[i],ProteinAlphabet()),id=proteinID[i],name=proteinID[i],description=description)
						fasta_out.append(out)
		self.sequences=fasta_out

	def writeFasta(self,output_file):
		fasta_out=self.sequences
		output_handle = open(output_file, "w")
		SeqIO.write(fasta_out, output_handle, "fasta")
		output_handle.close()


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("-gbk",help="gbk file, NCBI download",type=str,required=True)
	parser.add_argument("-out",help='Output filename',type=str,default='GBKextractor.fasta')
	args = parser.parse_args()
	GeneBank_file=args.gbk
	output_file=args.out

	fasta=ProteinSeqFromGenBank(GeneBank_file)
	fasta.writeFasta(output_file)
	print(str(len(fasta.sequences))+" Proteins extracted")

if __name__ == "__main__":
	main()