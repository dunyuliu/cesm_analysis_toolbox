for i in {1..16}
do 
  echo "Processing model $i ... ..."
  ncl case=$i enso_analysis.ncl 
done
