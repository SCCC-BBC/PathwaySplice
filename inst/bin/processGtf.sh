#!/bin/bash
tail -n +6 $1 | awk -F "\t" '{OFS="\t"; $1 = "chr"$1; print}' | awk -F"\t" '{OFS="\t"; if($1=="chrMT") $1="chrM"; print}' | sort -k1,1 -k4,4n > $2