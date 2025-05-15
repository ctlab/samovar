kraken2_samovar () {
  R1=($(ls -d *R1*))
  R2=($(ls -d *R2*))
  DB=~/kraken2_base
  
  mkdir $1
  for i in "${!R1[@]}"; do
    concat=$(echo ${R1[i]} | sed 's/_.*//g')
    
    echo -e "\n" $concat "\n"

    kraken2 \
      --use-names \
      --db $DB  \
      --threads 250 \
      --paired $2/${R1[i]} $2/${R2[i]}  \
      --report $1/$concat.report \
      --output $1/$concat.out 

  done

}

#kraken2_samovar reports/initial/k2_plus ./ /nfs/home/dsmutin/kraken2_bases
kraken2_samovar reports/initial/k2 ./