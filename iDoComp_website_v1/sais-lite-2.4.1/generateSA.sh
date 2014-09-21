#!/bin/bash
T="$(date +%s)"
rm -rf $1/*~
FILE=$1/*
ctr=0
for f in $FILE
do
  filename=$(basename "$f")
  extension="${filename##*.}"
  filename="${filename%.*}"
  
  echo "$f"
  tail -n +2 $f > $filename.tmp
  cat $filename.tmp | tr -d '\n' > tmp_chr$ctr.fa
  rm $filename.tmp

  ./sa.run tmp_chr$ctr.fa $2/$filename.sa > $2/$filename.time
 
  rm -rf tmp_chr$ctr.fa
 
  ctr=$[ctr + 1]

done
T="$(($(date +%s)-T))"
echo "Total time: ${T} seconds" > $2/total.time 
# ctr=$[ctr - 1]
# start=0
# 
# for i in $(eval echo "{$start..$ctr}")
# do
#   echo "$i"
#   ./STgen.run tmp_oneline_ref/chr$i.fa $2/$filename.tree > $2/$filename.time
# done

#rm tmp_oneline_tree/*
# rm -rf tmp_oneline_ref/*
