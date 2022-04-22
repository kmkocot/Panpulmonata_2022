#!/bin/bash
#SBATCH --job-name=TransDecoder
#SBATCH --mem=16G
#SBATCH -n 1 #tasks
#SBATCH -N 1 #nodes
#SBATCH -c 16 #cores per task
#SBATCH -o slurm_output-TransDecoder.%J
#SBATCH -e slurm_error-TransDecoder.%J
#SBATCH -p ultrahigh
#SBATCH --qos kmkocot
#SBATCH --mail-type=ALL
#SBATCH --mail-user=

export DK_ROOT=/share/apps/dotkit
. /share/apps/dotkit/bash/.dk_init
use hmmer3.1b2
use trinity
use transdecoder2.0.1
use bioinfoGCC
use perl5.22.1
use Anaconda

echo "Running TransDecoder.LongOrfs..."
for FILENAME in *.fa
do
TransDecoder.LongOrfs -t $FILENAME
done

echo "Running Diamond to search transcripts against uniprot/sprot..."
for FILENAME in ./*/longest_orfs.pep
do
DIR=`echo $FILENAME | cut -d "/" -f 1-2`
echo $FILENAME
echo $DIR
diamond blastp --threads 16 --query $FILENAME --db /kmk/databases/uniprot_sprot_21_Sept_2016/uniprot_sprot.fasta.dmnd --max-target-seqs 1 --evalue 0.00005 --outfmt 6 --out $DIR/longest_orfs.pep.blast_results
done

echo "Running TransDecoder.Predict..."
cp /kmk/databases/Pfam27/Pfam-A.hmm .
for FILENAME in *.fa
do
TransDecoder.Predict -t $FILENAME --retain_long_orfs 500 --retain_pfam_hits Pfam-A.hmm --retain_blastp_hits ./$FILENAME.transdecoder_dir/longest_orfs.pep.blast_results
done
rm Pfam-A.hmm



