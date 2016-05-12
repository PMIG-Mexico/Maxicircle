#!/bin/bash
cd AA/
for i in $(ls *.faa); 
do 
clustalo  --output-order=tree-order  -i $i -o $i.aln; 
makeblastdb -in $i -dbtype 'prot' -out $i ;
blastp -out $i.xml -outfmt 5 -query "/Users/Said/Github/Maxicircle/Leishmania_mexicana_maxicircle.fasta" -db $i ;
done;

mkdir ALN
mv *.aln ALN

cd ../DNA/

for i in $(ls *.fna); 
do 
clustalo  --output-order=tree-order  -i $i -o $i.aln; 
makeblastdb -in $i -dbtype 'nucl' -out $i ;
blastn -out $i.xml -outfmt 5 -query "/Users/Said/Github/Maxicircle/Leishmania_mexicana_maxicircle.fasta" -db $i ;
done;

mkdir ALN
mv *.aln ALN

cd ..
