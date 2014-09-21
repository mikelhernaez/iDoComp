#!/bin/bash

# This script generates the two .txt files needed for compression and decompression
# The input to the script are the following:

# $1: folder with the chromosomes of the reference genome
# $2: folder with the chromosomes of the target genome
# $3: folder with the trees of the reference genome
# $4: output folder for the reconstructed genomes
# $5: name of file for the encode.txt
# $6: name of file for the decode.txt

#### Important remark: ####
# The files in both the target and the reference folder should have the same order.

rm -f $5
rm -f $6

counterRef=0
for refFile in $1/* ; do
  counterRef=$[$counterRef +1]
  counterTar=0
  fileNameRef=$(basename "$refFile")
  fileNameRefExt="${fileNameRef##*.}"
  fileNameRef="${fileNameRef%.*}"

  # Write to decode.txt
  echo $refFile $4/${fileNameRef}_reconstructed.$fileNameRefExt >> $6

  for targetFile in $2/* ; do
    counterTar=$[$counterTar +1]
    if [ "$counterRef" = "$counterTar" ]; then

      # Write to encode.txt
      echo $refFile $targetFile $3/$fileNameRef.sa >> $5

    fi
  done
done

