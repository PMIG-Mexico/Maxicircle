cd AA/
for i in $(ls *.faa); 
do 
clustalo  --output-order=tree-order  -i $i -o $i.aln; 
done;

mkdir ALN
mv *.aln ALN

cd ../DNA/

for i in $(ls *.fna); 
do 
clustalo  --output-order=tree-order  -i $i -o $i.aln; 
done;

mkdir ALN
mv *.aln ALN

cd ..
