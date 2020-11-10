#!/bin/bash

val=0
while [ $val -le 10 ]
do
	if [ $val -eq 5 ];then
		((val=$val+2))
		continue;
	else
		echo "val=$val"
		((val++))
	fi
done
exit 0
