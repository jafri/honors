# Fetch data
prefetch --option-file acc.txt

# SRA to fastq
fastq-dump –X 5 –Z –split-files SRR000001