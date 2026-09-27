#!/bin/bash

df=$1
ids=$2

num=$(basename "$ids" .txt)

awk -F '\t' 'FNR == NR { ids[$1] = 1; next } $1 in ids { print $0 }' "$ids" "$df" > "$num.tsv"
