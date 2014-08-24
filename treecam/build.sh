#!/bin/bash -x

rm f*
./efficut/compressedcuts -r $1 -b16 -c1  -m1 -s8

rs_cnt=`ls -l f* | wc -l`

./treecam -n $rs_cnt

#files=`ls f*`
#num=0

#for file in $files
#do
#    ../hc/hypc -r $file -b16 -s8 -o "tree$num" -h2
#    num=$(($num+1))
#done
#mv -f tree* ../realrun/
#cd ../realrun/
#./realpc -n $num -r $1 -t $2 

