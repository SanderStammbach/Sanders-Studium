#!/bin/sh
for FILE in *.txt; do cat $FILE; done

ls -a ~/.ssh

for FILE in ~/.ssh/* ; do cat $FILE; done
