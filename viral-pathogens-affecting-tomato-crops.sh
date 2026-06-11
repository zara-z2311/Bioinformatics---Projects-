#Data retrieval
Source: NCBI Nucleotide Database
Server IP: 203.101.231.144
Directory path: /mnt/s-ws/everyone/projects/ToV
File name: ToCV_32_IDs.txt ToMV_69_IDs.txt
Project ID: ToV
cd /mnt/s-ws/everyone/projects/ToV/
source /mnt/s-ws/env/ToV_env.sh
cd
mkdir -p projects/Tov/
cd /mnt/s-ws/everyone/projects/ToV/
cp *txt ~/projects/Tov/
awk ‘NR%4==0 {print $1}’ ToMV_69_IDs.txt > ToMV_69_IDs.list
sed -n '4~4p' ToMV_69_IDs.txt |cut -d " " -f 1
efetch -db nuccore -input ToMV_69_IDs.list -format fasta > \
ToMV_69 _sequence.fasta

#Building Phylogenetic Tree

mafft --localpair --maxiterate 1000 --adjustdirection ToMV_69_sequence.fasta > ToMV_69_mafft.fasta
trimal -in ToMV_69_mafft.fasta -nexus -out ToMV_69_mafft.nex

#Generate MrBayes Tree

mb
begin mrbayes;
set autoclose=no nowarn=yes;
lset nst=6 rates=invgamma;
mcmcp ngen= 5000000 printfreq=1000 samplefreq=100 nchains=4 savebrlens=yes;
mcmc;
end;
exe ToMV_69_mafft.nex
mcmc
sump
sumt

#Translate genome fasta into proteins

transeq ToCV_32_combined.fasta ToCV_32_seq.prot.fas
getorf ToCV_32_combined.fasta ToCV_32_orf.fasta -minsize 400
cat ToCV_32_orf.fasta | grep "^>" | cut -d " " -f 1 | sed 's/_[0-9]*$//g' | sort | uniq -c
cat ToMV_69_orf.fasta | grep "^>" | cut -d " " -f 1 | sed 's/_[0-9]*$//g' | sort | uniq -c
efetch -db nuccore -input ToCV_32_IDs.list -format fasta_cds_aa > ToCV_proteins.fasta
grep ">" ToCV_proteins.fasta | head -20
cat ToMV_69_orf.fasta | grep "^>" | cut -d " " -f 1 | sed 's/_[0-9]*$//g' | sort | uniq -c > ToMV_orf_counts.txt
column -t ToMV_orf_counts.txt | head
awk '{print $1}' ToMV_orf_counts.txt | sort -n | uniq -c
grep ">" ToCV_32_orf.fasta | cut -d"_" -f1 | sort | uniq
getorf ToCV_32_sequence.fasta ToCV_32_orf_min200.fasta -minsize 200
cat ToCV_32_orf_min200.fasta | grep "^>" | cut -d" " -f1 | sed 's/_[0-9]*$//g' | sort | uniq -c
grep ">" ToCV_proteins.fasta | head -20
cat ToCV_32_orf_min200.fasta | grep "^>" | cut -d " " -f 1 | sed 's/_[0-9]*$//g' | sort | uniq -c > ToCV_orf_min200_counts.txt
column -t ToCV_orf_min200_counts.txt
awk '{print $1}' ToCV_orf_min200_counts.txt | sort -n | uniq -c

#Identify the top conserved proteins and find known protein fucntions

makeblastdb -in ToCV_32_seq.prot.fas -dbtype prot -out ToCV_db
blastp -query ToMV_seq.prot.fast -db ToCV_db -out ToMV_CV_blast.txt -outfmt 6 -evalue 1e-1 -max_target_seqs 10
seqkit grep -f ToMV_top.list ToMV_proteins.fasta |seqkit subseq -r 1175:1544 > ToMV_list.fasta
cat ToMV_CV_blast.txt|cut -d "_" -f 1|sort|uniq -c > ToMV_conserved.txt
