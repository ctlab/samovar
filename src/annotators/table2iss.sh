path=/data/test_genomes

mkdir tests_outs/benchmarking/initial
cd tests_outs/benchmarking/initial

iss_generate () {
  iss generate \
    --genomes $1 \
    --model hiseq \
    --output ${2}_cont \
    --n_reads $3
}

iss_generate_all () {
    for i in $(ls $path/*.fna); do
        iss_generate $i $(basename $i .fna)
    done
}

for ((i=1;i<=25;i++)); do
  # making network
  human_reads=$((($RANDOM % 100) * 250))
  meta_reads=$(((25000 - $human_reads)/3))
  
  # generate
  iss_generate_all $i
  
  # concat
  cat *${i}*_R1* > ${i}_full_R1.fastq
  cat *${i}*_R2* > ${i}_full_R2.fastq

  # clean
  rm *tmp*
  rm *_cont_*
  rm *abundance*
  rm *vcf
  
  echo -e "done " $i "\n"
done