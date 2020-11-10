#!/bin/bash

cd /home/test/linux/
unzip bash_homework.zip
cd bash_homework/
dir=$(ls)
for val in $dir
do
	if [ -f $val ];then
		echo $val >> ../filenames.txt
	elif [ -d $val ];then
		echo $val >> ../dirnames.txt
	fi
done
exit 0
