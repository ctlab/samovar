kaiju_samovar () {
  DB=~/kaiju_refseq

  R1=($(ls -d *R1*))
  R2=($(ls -d *R2*))
  
  mkdir $1
  
  for i in "${!R1[@]}"; do
    concat=$(echo ${R1[i]} | sed 's/_.*//g')
    
    echo -e "\n" $concat "\n"

    kaiju \
      -t $DB/nodes.dmp \
      -f $DB/kaiju_db_refseq_ref.fmi \
      -i $2/${R1[i]} \
      -j $2/${R2[i]} \
      -z 25000 \
      -o $1/$concat.raw
      
    kaiju-addTaxonNames \
    -t $DB/nodes.dmp \
    -n $DB/names.dmp \
    -i $1/$concat.raw \
    -o $1/$concat.out
    
  done
}

kaiju_samovar reports/initial/kaiju ./