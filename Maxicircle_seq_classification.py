# -*- coding: utf-8 -*-
"""
Created on Thu May  11 2016

@author: Said
"""
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import ProteinAlphabet
import argparse
import os

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
						out=SeqRecord(seq=Seq(proteinSequence[i],ProteinAlphabet()),id=proteinID[i],name=product[i],description=description)
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

def writeDNA(fasta_out,output_file):
	output_file=output_file+'.genes.fna'
	output_handle = open(output_file, "w")
	SeqIO.write(fasta_out, output_handle, "fasta")
	output_handle.close()

def writeAA(fasta_out,output_file):
	output_file=output_file+'.proteins.faa'
	output_handle = open(output_file, "w")
	SeqIO.write(fasta_out, output_handle, "fasta")
	output_handle.close()


class sequenceClassification(FeatureFromGenBank):
	def __init__(self,seq,out_dir):
		genes=seq.sequencesDNA
		genesID=[x.id for x in genes]
		uniquegenes=list(set([x.upper() for x in genesID]))
		current_path=os.getcwd()
		try:
			os.mkdir(out_dir)
			os.mkdir(out_dir+"/DNA")
		except OSError as exc:
			pass
		os.chdir(out_dir+"/DNA/")
		for ID in uniquegenes:
			indices = [i for i, x in enumerate(genesID) if x == ID]
			genSeq=[]
			genSeq=[genes[int(j)] for j in indices]
			writeDNA(genSeq, ID.replace(" ","_"))
		os.chdir(current_path)
		proteins=seq.sequencesAA
		protID=[x.name for x in proteins]
		uniqueprot=list(set([x.upper() for x in protID]))
		try:
			os.mkdir(out_dir+"/AA")
		except OSError as exc:
			pass		
		os.chdir(out_dir+"/AA/")
		for ID in uniqueprot:
			indices = [i for i, x in enumerate(protID) if x == ID]
			protSeq=[]
			protSeq=[proteins[int(j)] for j in indices]
			writeAA(protSeq, ID.replace(" ","_"))
		os.chdir(current_path)
		print(str(len(uniqueprot))+" Differents proteins classified \n")
		print(str(len(uniquegenes))+" Differents genes classified \n")




def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("-gbk",help="gbk file, NCBI download",type=str,required=True)
	parser.add_argument("-out_dir",help='Output filename',type=str,default='Classification')
	args = parser.parse_args()
	GeneBank_file=args.gbk
	out_dir=args.out_dir

	fasta=FeatureFromGenBank(GeneBank_file)
	print(str(len(fasta.sequencesAA))+" Proteins extracted\n")
	print(str(len(fasta.sequencesDNA))+" Genes extracted")
	sequenceClassification(fasta, out_dir)

if __name__ == "__main__":
	main()