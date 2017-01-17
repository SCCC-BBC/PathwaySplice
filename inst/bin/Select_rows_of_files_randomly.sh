#!/bin/bash

input_file=$1
num_lines=$2
output_file=$3

gshuf -n "$num_lines" "$input_file" > "$output_file" 