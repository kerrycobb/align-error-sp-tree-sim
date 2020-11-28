
for i in {1..200}; do
  DIR=out-sp2-gen4-loc10000-len1000/seed$i-reps1/filter-0.01-upper/
  echo $DIR
  ls $DIR | wc -l
  # WC=$(ls $DIR | wc -l)
  # echo $WC

  # END=$(tail -1 $i)
  # if [[ $END != "Simulation Complete" ]]; then
  #   echo $i
    #SCRIPT=$(cut -d' ' -f11- <<<$END)
    #rm $i
    #myqsub -N sim-250-$i -t 2:00:00 -m 10gb \"$SCRIPT\"
done
