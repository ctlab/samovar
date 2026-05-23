metaphlan_samovar () {
  R1=($(ls -d [1-9]_*R1*))
  R2=($(ls -d *R2*))
  
  mkdir $1
  
  for i in "${!R1[@]}"; do
    concat=$(echo ${R1[i]} | sed 's/_.*//g')
    
    echo -e "\n" $concat "\n"

    metaphlan \
      --input_type fastq \
      --nproc 2500 \
      $2/${R1[i]} \
      --bowtie2out metagenome.bowtie2.bz2 \
      -o $1/$concat.out 
      
      rm metagenome*
  done
}

metaphlan_samovar reports/initial/metaphlan ./