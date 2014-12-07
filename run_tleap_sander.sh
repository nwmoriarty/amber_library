#! /bin/bash


for i in 0 1 2 3 4 5 6 7 8 9 a b c d e f g h i j k l m n o p q r s t u v w x y z
do
  cd $i

  n=`ls -1 *mol2 |wc -l`
  count=0
  for i in `ls -1 *.mol2`; do
    let count=count+1
#    if [ $count -eq 3 ]; then cd ..; exit; fi
    code=${i%.mol2}
    echo $count of $n:   $code
    ../tleap_sander.sh $code
    ../tleap_sander_igb.sh $code
  done
  cd ..
done


  
