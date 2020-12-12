#!/bin/bash

# this script will generate a file with orthogroups that have duplicates in both V. appalachiana and V. lineata and contains one copy of the outgroup A. platyneuron 

for FILE in *.fa; do 
	NDUV=$(grep -o "NDUV" $FILE | wc -w)
	SKYV=$(grep -o "SKYV" $FILE | wc -w)
	KJZG=$(grep -o "KJZG" $FILE | wc -w)
	if [ "$NDUV" != 4 ]; then NDUV=NOTDUP; fi #there are 2 taxa codes per header
	if [ "$SKYV" != 4 ]; then SKYV=NOTDUP; fi #there are 2 taxa codes per header
	if [ "$KJZG" != 2 ]; then KJZG=NOTDUP; fi #there are 2 taxa codes per header
	all=$(echo $NDUV $SKYV $KJZG)
	if [[ $all == *"NOTDUP"* ]]; then echo "No duplicates, end!"; else echo "$FILE" >> vittaria_duplicates; fi;
done
