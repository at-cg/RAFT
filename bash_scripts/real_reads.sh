wget -i README
cat *.fastq.1 > reads.fastq
rm *.fastq.1
sed -n '1~4s/^@/>/p;2~4p' reads.fastq > reads.fasta
python $SCRATCH/code/add_read_ids_fasta.py reads.fasta output_reads.fasta
rm reads.fasta
rm reads.fastq