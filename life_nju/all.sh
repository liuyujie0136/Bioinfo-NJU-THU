#!/bin/bash

if [ all.out ];then rm all.out; fi

dir=`ls -F | grep "/$"`
for i in $dir
do
	echo $i >> all.out
	cd $i
	files=`ls`
	for j in $files
	do
		dos2unix $j
		echo $j >> ../all.out
		awk 'BEGIN{ORS=" "}; /^[^0-9]/{print}' $j >> ../all.out
		echo -e "\n" >> ../all.out
	done
	cd ..
done
exit 0
