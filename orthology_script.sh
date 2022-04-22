#!/bin/bash
#
# Kevin M. Kocot
# 
#This script takes the output of OrthoFinder and performs several "cleaning" steps to remove groups and sequences that are not suitable for phylogenomic analysis.
#Note which programs must be included in $PATH
#Run this script from the working directory containing the input fasta files (extension must be .fa) where each file represents one putatively orthologous group.
#Fasta headers must be in the following format: >species_name_abbreviation|annotation_or_sequence_ID_information
#Example: >LGIG|Contig1234
#Fasta headers may not include spaces or non-alphanumeric characters except for underscores (pipes are OK as field delimiters only).


########################################################################
################## Change these before your first use ##################
########################################################################
#Set values for variables.
MIN_SEQUENCE_LENGTH=100 #Deletes original seqs shorter than this length
MIN_ALIGNMENT_LENGTH=100 #Minimum length of a trimmed alignment in amino acids
HALF_TAXA=45 #Minimum number of OTUs to keep an OG
THREEFOURTHS_TAXA=68
CLASS_PATH=/usr/local/bin/ #Location of PhyloTreePruner .class files
CORES=35 #Number of CPU cores to use with Parallel
PREFIX=IQ-TREE2 #IQ-TREE 2 prefix for output files
########################################################################


#Fix headers and remove newlines.
echo "Removing linebreaks in sequences and fixing headers..."
for FILENAME in *.fa
do
sed -i ':a; $!N; /^>/!s/\n\([^>]\)/\1/; ta; P; D' $FILENAME
done
echo Done
echo


#Backup all sequences
echo "Making a backup of all sequences before beginning..."
mkdir unedited_sequences
cp *0.fa ./unedited_sequences/
cp *1.fa ./unedited_sequences/
cp *2.fa ./unedited_sequences/
cp *3.fa ./unedited_sequences/
cp *4.fa ./unedited_sequences/
cp *5.fa ./unedited_sequences/
cp *6.fa ./unedited_sequences/
cp *7.fa ./unedited_sequences/
cp *8.fa ./unedited_sequences/
cp *9.fa ./unedited_sequences/
echo
echo Done
echo


#Delete sequences shorter than $MIN_SEQUENCE_LENGTH
echo "Deleting sequences shorter than $MIN_SEQUENCE_LENGTH AAs..."
for FILENAME in *.fa
do
grep -B 1 "[^>].\{$MIN_SEQUENCE_LENGTH,\}" $FILENAME > $FILENAME.out
sed -i 's/--//g' $FILENAME.out
sed -i '/^$/d' $FILENAME.out
rm -rf $FILENAME
mv $FILENAME.out $FILENAME
done
echo Done
echo


#If fewer than $HALF_TAXA different species are represented in the file, move that file to the "rejected_few_taxa" directory. 
echo "Removing groups with fewer than $HALF_TAXA taxa..."
mkdir -p rejected_few_taxa_1
for FILENAME in *.fa
do
awk -F"|" '/^>/{ taxon[$1]++ } END{for(o in taxon){print o,taxon[o]}}' $FILENAME > $FILENAME\.taxon_count #Creates temporary file with taxon abbreviation and number of sequences for that taxon in $FILENAME
taxon_count=`grep -v 0 $FILENAME\.taxon_count | wc -l` #Counts the number of lines with an integer >0 (= the number of taxa with at least 1 sequence)
if [ "$taxon_count" -lt "$HALF_TAXA" ] ; then 
echo $FILENAME
mv $FILENAME ./rejected_few_taxa_1/
fi
done
rm -rf *0.fa.taxon_count
rm -rf *1.fa.taxon_count
rm -rf *2.fa.taxon_count
rm -rf *3.fa.taxon_count
rm -rf *4.fa.taxon_count
rm -rf *5.fa.taxon_count
rm -rf *6.fa.taxon_count
rm -rf *7.fa.taxon_count
rm -rf *8.fa.taxon_count
rm -rf *9.fa.taxon_count
echo Done
echo


#Remove redundant sequences using uniqHaplo
mkdir preUniqHaplo
cp *.fa preUniqHaplo
echo "Removing redundant sequences using uniqHaplo..."
ls *0.fa | parallel -j $CORES 'perl /usr/bin/uniqHaplo.pl -a {} > {}.uniq'
ls *1.fa | parallel -j $CORES 'perl /usr/bin/uniqHaplo.pl -a {} > {}.uniq'
ls *2.fa | parallel -j $CORES 'perl /usr/bin/uniqHaplo.pl -a {} > {}.uniq'
ls *3.fa | parallel -j $CORES 'perl /usr/bin/uniqHaplo.pl -a {} > {}.uniq'
ls *4.fa | parallel -j $CORES 'perl /usr/bin/uniqHaplo.pl -a {} > {}.uniq'
ls *5.fa | parallel -j $CORES 'perl /usr/bin/uniqHaplo.pl -a {} > {}.uniq'
ls *6.fa | parallel -j $CORES 'perl /usr/bin/uniqHaplo.pl -a {} > {}.uniq'
ls *7.fa | parallel -j $CORES 'perl /usr/bin/uniqHaplo.pl -a {} > {}.uniq'
ls *8.fa | parallel -j $CORES 'perl /usr/bin/uniqHaplo.pl -a {} > {}.uniq'
ls *9.fa | parallel -j $CORES 'perl /usr/bin/uniqHaplo.pl -a {} > {}.uniq'
rm -rf *.fa
rename 's/.fa.uniq/.fa/g' *.fa.uniq
echo Done
echo


#Align the remaining sequences using Mafft.
echo "Aligning sequences using Mafft (auto)..."
mkdir backup_alignments
ls *0.fa | parallel -j $CORES 'mafft --auto --reorder --quiet --localpair --maxiterate 1000 {} > {}.aln'
ls *1.fa | parallel -j $CORES 'mafft --auto --reorder --quiet --localpair --maxiterate 1000 {} > {}.aln'
ls *1.fa | parallel -j $CORES 'mafft --auto --reorder --quiet --localpair --maxiterate 1000 {} > {}.aln'
ls *2.fa | parallel -j $CORES 'mafft --auto --reorder --quiet --localpair --maxiterate 1000 {} > {}.aln'
ls *3.fa | parallel -j $CORES 'mafft --auto --reorder --quiet --localpair --maxiterate 1000 {} > {}.aln'
ls *4.fa | parallel -j $CORES 'mafft --auto --reorder --quiet --localpair --maxiterate 1000 {} > {}.aln'
ls *5.fa | parallel -j $CORES 'mafft --auto --reorder --quiet --localpair --maxiterate 1000 {} > {}.aln'
ls *6.fa | parallel -j $CORES 'mafft --auto --reorder --quiet --localpair --maxiterate 1000 {} > {}.aln'
ls *7.fa | parallel -j $CORES 'mafft --auto --reorder --quiet --localpair --maxiterate 1000 {} > {}.aln'
ls *8.fa | parallel -j $CORES 'mafft --auto --reorder --quiet --localpair --maxiterate 1000 {} > {}.aln'
ls *9.fa | parallel -j $CORES 'mafft --auto --reorder --quiet --localpair --maxiterate 1000 {} > {}.aln'
rm -rf *0.fa
rm -rf *1.fa
rm -rf *2.fa
rm -rf *3.fa
rm -rf *4.fa
rm -rf *5.fa
rm -rf *6.fa
rm -rf *7.fa
rm -rf *8.fa
rm -rf *9.fa
rename 's/.fa.aln/.fa/g' *.fa.aln
cp *.fa ./backup_alignments/
echo Done
echo


#Remove newlines.
echo "Removing linebreaks in sequences..."
for FILENAME in *.fa
do
sed -i ':a; $!N; /^>/!s/\n\([^>]\)/\1/; ta; P; D' $FILENAME
done
echo Done
echo


#Clean alignments with HmmCleaner
echo "Removing misaligned sequence regions with HmmCleaner..."
mkdir HmmCleaner_files
ls *.fa | parallel -j $CORES 'HmmCleaner.pl {} --specificity'
mv *.fa ./HmmCleaner_files
mv *.log ./HmmCleaner_files
mv *.score ./HmmCleaner_files
rename 's/_hmm.fasta/.fa/g' *_hmm.fasta
echo Done
echo


#Trim alignments with BMGE
echo "Trimming ambiguously aligned columns in the alignment with BMGE..."
mkdir backup_pre-BMGE
ls *.fa | parallel -j $CORES 'java -jar /usr/bin/BMGE-1.12/BMGE.jar -i {} -t AA -of {}.BMGE'
mv *.fa ./backup_pre-BMGE
cp *.BMGE ./backup_pre-BMGE
rename 's/.BMGE//g' *.BMGE
echo Done
echo


#Remove newlines.
echo "Removing linebreaks in sequences..."
for FILENAME in *.fa
do
sed -i ':a; $!N; /^>/!s/\n\([^>]\)/\1/; ta; P; D' $FILENAME
done
echo Done
echo


#Remove any sequences that don't overlap with all other sequences by at least 20 amino acids. This check runs until all sequences overlap with all other sequences by at least 20 amino acids.
echo "Removing short sequences that don't overalp with all other sequences by at least 20 AAs..."
for FILENAME in *.fa
do
java -cp /usr/bin AlignmentCompare $FILENAME
done
echo Done
echo
rm -rf myTempFile.txt


#If fewer than $HALF_TAXA different species are represented in the file, move that file to the "rejected_few_taxa" directory. 
echo "Removing groups with fewer than $HALF_TAXA taxa..."
mkdir -p rejected_few_taxa_2
for FILENAME in *.fa
do
awk -F"|" '/^>/{ taxon[$1]++ } END{for(o in taxon){print o,taxon[o]}}' $FILENAME > $FILENAME\.taxon_count #Creates temporary file with taxon abbreviation and number of sequences for that taxon in $FILENAME
taxon_count=`grep -v 0 $FILENAME\.taxon_count | wc -l` #Counts the number of lines with an integer >0 (= the number of taxa with at least 1 sequence)
if [ "$taxon_count" -lt "$HALF_TAXA" ] ; then 
echo $FILENAME
mv $FILENAME ./rejected_few_taxa_2/
fi
done
rm -rf *.fa.taxon_count
echo Done
echo


#Makes a tree for each OG with FastTree
echo "Making a tree for each OG using FastTreeMP..."
for FILENAME in *.fa
do
FastTreeMP -slow -gamma $FILENAME > $FILENAME.tre
done
rename 's/.fa.tre/.tre/g' *.fa.tre
echo Done


#Runs PhyloPyPruner and makes trees
echo "Running PhyloPyPruner, FastTree 2, and IQ-TREE 2..."

phylopypruner --threads $CORES --dir . --min-taxa $HALF_TAXA --min-support 0.9 --mask pdist --trim-lb 3 --trim-divergent 0.75 --min-pdist 0.01 --prune LS #Run with version 0.9.5
cd phylopypruner_output
FastTreeMP -slow -gamma supermatrix.fas > FastTree.tre
iqtree2 -T AUTO -B 1000 -s supermatrix.fas -Q partition_data.txt --prefix $PREFIX -m MFP --wbtl
cd ..
mv phylopypruner_output phylopypruner_output_min_50_pct

phylopypruner --threads $CORES --dir . --min-taxa $THREEFOURTHS_TAXA --min-support 0.9 --mask pdist --trim-lb 3 --trim-divergent 0.75 --min-pdist 0.01 --prune LS #Run with version 0.9.5
cd phylopypruner_output
FastTreeMP -slow -gamma supermatrix.fas > FastTree.tre
iqtree2 -T AUTO -B 1000 -s supermatrix.fas -Q partition_data.txt --prefix $PREFIX -m MFP --wbtl
cd ..
mv phylopypruner_output phylopypruner_output_min_75_pct

echo Done





