#!/bin/bash

declare -i sum
sum=0
for i in `seq 1 2 100`
do
	sum=$sum+$i
done
echo $sum
exit 0
