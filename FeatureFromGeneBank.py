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

class FeatureFromGenBank:
	def __init__(self,GeneBank_file):
		fasta_out=[]
		genes_out=[]
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
						product=['unknown']*len(proteinSequence)
					for i in range(len(proteinSequence)):
						description='|'+organism+'|'+product[i]
						out=SeqRecord(seq=Seq(proteinSequence[i],ProteinAlphabet()),id=proteinID[i],name=proteinID[i],description=description)
						fasta_out.append(out)
				if 'gene' in feature.type:
					startgene=int(feature.location.start)
					endgene=int(feature.location.end)
					if 'gene' in feature.qualifiers:
						geneID="".join(feature.qualifiers['gene'])
					elif 'db_xref' in feature.qualifiers:
						geneID=" ".join(feature.qualifiers['db_xref'])
					else:
						geneID='unknown'
					description='|'+organism+'|'+name
					gout=SeqRecord(seq=seq_record.seq[startgene:endgene],id=geneID,name=geneID,description=description)
					genes_out.append(gout)
		self.sequencesAA=fasta_out
		self.sequencesDNA=genes_out


	def writeFAA(self,output_file):
		fasta_out=self.sequencesAA
		output_file=output_file+'.proteins.faa'
		output_handle = open(output_file, "w")
		SeqIO.write(fasta_out, output_handle, "fasta")
		output_handle.close()

	def writeDNA(self,output_file):
		fasta_out=self.sequencesDNA
		output_file=output_file+'.genes.fna'
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

	fasta=FeatureFromGenBank(GeneBank_file)
	fasta.writeFAA(output_file)
	fasta.writeDNA(output_file)
	print(str(len(fasta.sequencesAA))+" Proteins extracted\n")
	print(str(len(fasta.sequencesDNA))+" Genes extracted")

if __name__ == "__main__":
	main()