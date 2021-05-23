#Manu Singh 26 Dec 2021

export PATH=/programs/salmon-1.2.1/bin:$PATH

my_index = /workdir/Manu/Homo_Genome/Genome/bed_files/Alevin/Genes_TEindex

my_tgmap = /workdir/Manu/Homo_Genome/Genome/bed_files/Alevin/TE_long/txp2gene.tsv

# Running Alevin on the transcriptome that we build using genes and TEs

for i in $(ls ./10XDISCEL*  | sed s/_[12].fastq.gz// | sort -u)
do /programs/salmon-1.2.1/bin/salmon alevin -l ISR \
-1 ${i}_1.fastq.gz \
-2 ${i}_2.fastq.gz  \
--chromium  \
-i $my_index \
-p 24 \
-o quants/${i}_quant \
--tgMap $my_tgmap \
--expectCells 10000 \
--forceCells 10000
 done

cd quants

for old in *
        do 
           if [ -f "$old" ]
              then new=$(echo $old | sed 's/10XDISCEL//g'); mv "$old" "$new"
           fi
        done



#Each output for each sample was then renamed manually with their original samples identities