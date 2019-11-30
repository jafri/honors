#! /usr/bin/env bash
while read p; do
  echo "#! /usr/bin/env bash
  #BSUB -J fdmp_${p}
  #BSUB -n 1
  #BSUB -R "span[hosts=1]"
  #BSUB -R "rusage[mem=3500]"
  #BSUB -R "select[type==any]"
  #BSUB -M 4000
  #BSUB -We 00:10
  #BSUB -W 02:00
  #BSUB -o test_%J.out
  #BSUB -e test_%J.err

  kallisto quant --index=/gpfs/home/sjafri/TEST/kallisto_index --output-dir=/gpfs/home/sjafri/scpTest/honors/${p} --threads=4 /gpfs/home/sjafri/scpTest/honors/${p}_1.fastq.gz /gpfs/home/sjafri/scpTest/honors/${p}_2.fastq.gz

  sleep 10
  echo 'finished with exit code $? at: `date`'" | bsub
done <acc.txt