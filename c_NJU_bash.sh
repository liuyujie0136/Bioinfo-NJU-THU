#!/bin/bash

cd /home/test/share/c_NJU/
touch ../c_NJU.md
files=`ls`
for i in $files
do
	echo '##' $i '\n```c++' >> ../c_NJU.md
	cat $i >> ../c_NJU.md
	echo '```\n' >> ../c_NJU.md
done
exit 0