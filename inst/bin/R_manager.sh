#!/bin/bash

#example:

#sh R_manager.sh track R/
 
#sh R_manager.sh exclude R/ GotermAnalysis2GetAllGOTerms_AdjustedByNumOfJunctionWithinOneGene.R

way=$1
DIR=$2
pattern=$3

echo 1.$way
echo 2.$DIR
echo 3.$pattern

for file in $(ls $DIR$pattern)

do

f=`echo "$file"`

if [ $way == 'exclude' ]
then

echo "exclude file"

echo "^"$f >> .Rbuildignore 

fi

done
