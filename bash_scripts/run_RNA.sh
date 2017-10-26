#!/bin/bash -e

#SBATCH --cpus-per-task=4
#SBATCH --mem=50G
#SBATCH --array=1-50%5

module load bedtools
module load bedops

sample=$( sed -n ${SLURM_ARRAY_TASK_ID}p samples.txt )

srun wget http://egg2.wustl.edu/roadmap/data/byDataType/rna/signal/unnormalized_wig/stranded/$sample\.gz
srun gunzip $sample\.gz
srun wig2bed < $sample > $sample\.bed
rm $sample
for j in refseq_exons.txt.sorted rmsk_TEother.txt.sorted;
do
    TEs=$j
    suffix=${TEs%.txt}
    srun bedtools intersect -wo -a $TEs -b $sample\.bed > $suffix\_$sample\.bed -sorted
    srun python /scratch/twlab/dli/ecp/calculate_average_RNA.py $suffix\_$sample\.bed $TEs $suffix\_$sample\_average.txt
    rm $suffix\_$sample\.bed
done
rm $sample\.bed
