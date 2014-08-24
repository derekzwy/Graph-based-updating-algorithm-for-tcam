#!/bin/bash -x

rm -rf f*
rm -rf tree*
./compressedcuts -r $1 -b16 -c1  -m1 -s8
files=`ls f*`
num=0

#for file in $files
#do
#    ../hc/hypc -r $file -b16 -s8 -o "tree$num" -h2
#    num=$(($num+1))
#done
#mv -f tree* ../realrun/
#cd ../realrun/
#./realpc -n $num -r $1 -t $2 

