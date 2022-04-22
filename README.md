# snpad-nf
Nextflow pipeline to run snpAD. 

```bash

export TAB=$(printf "\t")
cat > input.tsv <<EOF
Sample_Name${TAB}RG${TAB}BAM
Sample${TAB}ILLUMINA-Sample_L001,ILLUMINA-Sample_L002${TAB}/path/to/Sample.sorted.bam
EOF

# Assumes that snpad container image is extracted at /opt/images-live/${USER}
# and you are running charliecloud containers and the slurm scheduler

nextflow run \
	sc13-bioinf/snpad-nf \
	-profile charliecloud,slurm \
	--input input.tsv \
	--reference /opt/reference/human/chm13/chm13v2.0_maskedY_rCRS.fa \
	--chromosomes chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM \
	--prior_chromosomes chr21,chrX,chrY,chrM

```


