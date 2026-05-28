nextflow run nf-core/rnaseq \
  -r 3.9 \
  -profile icbi,singularity \
  -c rnaseq.config \
  -w /data/scratch/kvalem/projects/2021/microbial-metabolites/work
  --gencode
