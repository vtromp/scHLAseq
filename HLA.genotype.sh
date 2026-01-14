#
zcat $GTF | \
awk -F '\t' 'BEGIN{OFS="\t"}
$1=="chr6" &&
$3=="gene" &&
$9 ~ /gene_name "HLA-(DP|DQ|DR)/ {
  match($9, /gene_name "([^"]+)"/, g)
  print $1, $4-1, $5, g[1], ".", $7
}' > HLA_classII_genes.bed

#
samtools view -h -L HLA_classII_genes.bed $BAM | awk '$0 ~ /^@/ || ($0 ~ /\tCB:Z:/ && $0 ~ /\tUB:Z:/)' > HLA_classII_mapped_reads.sam
samtools view -h -f 4 $BAM | awk '$0 ~ /^@/ || ($0 ~ /\tCB:Z:/ && $0 ~ /\tUB:Z:/)' > unmapped_reads.sam
#
samtools view -b -o HLA_classII_mapped_reads.bam HLA_classII_mapped_reads.sam
samtools view -b -o unmapped_reads.bam unmapped_reads.sam
#
samtools merge -o HLA_classII_and_unmapped_reads.bam HLA_classII_mapped_reads.bam unmapped_reads.bam
#
rm HLA_classII_mapped_reads.sam unmapped_reads.sam HLA_classII_mapped_reads.bam unmapped_reads.bam

#
samtools fastq HLA_classII_and_unmapped_reads.sam > HLA_classII_and_unmapped_reads.fastq

#
samtools view HLA_classII_and_unmapped_reads.sam | \
awk '{
  read_id=$1;
  chr=$3;
  pos=$4;
  cb="NA";
  ub="NA";
  for(i=12;i<=NF;i++){
    if($i ~ /^CB:Z:/) { cb=substr($i,6) }
    if($i ~ /^UB:Z:/) { ub=substr($i,6) }
  }
  if(cb != "NA" && ub != "NA") {
    print read_id"\t"cb"\t"ub"\t"chr"\t"pos
  }
}' > overview_reads.tsv

#
prepare.HLA.database.R --class-II-only

#
vsearch --usearch_global HLA_classII_and_unmapped_reads.fastq --db HLA_alleles.fasta --id 0.96 --query_cov 1.0 --userfields query+target+id+mism+gaps+tilo+tihi --userout all_HLA_alignments.tsv

