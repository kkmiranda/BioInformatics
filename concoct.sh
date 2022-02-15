###Concoct
module load Anaconda3/2019.03
source activate concoct

for SET in MG2 MG4 MG5 MG6 MG7 MG8
do

#FASTA -> BEDfile (with cut_up_fasta.py into chunks of 10kb)
#In order to give more weight to larger contigs and mitigate the effect of assembly errors (LENGTH = 1000)
cut_up_fasta.py 02_FASTA/${SET}/${SET}-contigs-prefix-formatted-only.fa -c 10000 -o 0 --merge_last -b concoct_output/${SET}/${SET}_contigs_10K.bed > concoct_output/${SET}/${SET}_contigs_10K.fa

#BEDfile + BAMfile -> Coverage_file.tsv (with concoct_coverage_table.py)
concoct_coverage_table.py concoct_output/${SET}/${SET}_contigs_10K.bed 04_MAPPING/${SET}/${SET}.bam > concoct_output/${SET}/${SET}_coverage_table.tsv

#FASTA + coverage_file.tsv -> Autobinning with Concoct
concoct --coverage_file concoct_output/${SET}/${SET}_coverage_table.tsv --composition_file concoct_output/${SET}/${SET}_contigs_10K.fa -c 100 -l 1000 -b concoct_output/${SET}/

#Merge
merge_cutup_clustering.py concoct_output/${SET}/clustering_gt1000.csv > concoct_output/${SET}/${SET}_concoct.csv

##Reformat csv file to tab delimited text file
# replaces commas with tabs (to put in a tab in bash, type ctrl+v and then hit tab)
sed 's/,/	/g' concoct_output/${SET}/${SET}_concoct.csv > concoct_output/${SET}/${SET}_concoct.txt

#delete headers
sed -i 1d concoct_output/${SET}/${SET}_concoct.txt 

#starts the second column with Bin_
awk '$2="Bin_"$2' concoct_output/${SET}/${SET}_concoct.txt > concoct_output/${SET}/${SET}_concoct_wBin.txt

#Converts it completely into tab delimited text
awk '{ for(i=1;i<=NF;i++){if(i==NF){printf("%s\n",$NF);}else {printf("%s\t",$i)}}}' concoct_output/${SET}/${SET}_concoct_wBin.txt > concoct_output/${SET}/${SET}_concoct.txt

module load Anaconda3/2019.03
conda activate anvio-6.1

anvi-import-collection concoct_output/${SET}/${SET}_concoct.txt \
                       -c 03_CONTIGS/${SET}-contigs.db \
                       -p 05_ANVIO_PROFILE/${SET}/${SET}/PROFILE.db \
                       -C ${SET}_CONCOCT \
                       --contigs-mode

anvi-show-collections-and-bins -p 05_ANVIO_PROFILE/${SET}/${SET}/PROFILE.db

done

# anvi-delete-collection -p 05_ANVIO_PROFILE/MG3/MG3/PROFILE.db -C MG3_CONCOCT #In case this is necessary.
# anvi-interactive -c 03_CONTIGS/MG3-contigs.db -p 05_ANVIO_PROFILE/MG3/MG3/PROFILE.db -C MG3_CONCOCT --server-only -P 8080

