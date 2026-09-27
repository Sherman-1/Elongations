#!/bin/bash
# ssearch.sh - Proteic alignment wrapper for ssearch36
# Tools paths can be set via environment variables or will use PATH defaults

query=$1
subjects=$2
ncpus=$3
mem=$4

# Use environment variable if set, otherwise use tool name (expects it in PATH)
SSEARCH=${SSEARCH:-ssearch36}

echo -e "qseqid\tqlen\tsseqid\tslen\tppos\talign_length\tmismatches\tgaps\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tCIGAR" > temp
$SSEARCH -3 -p -s BL50 -f -11 -g -1 -T${ncpus} -XM${mem}G -m8BCL -m "FA ${query}_vs_neighbors.aln" $query $subjects >> temp

/usr/bin/python3 -c "
import polars as pl
(
	pl.scan_csv('temp', separator = '\t', has_header = True)
	.with_columns(
		qcovhsp = ((pl.col('align_length') / pl.col('qlen')) * 100)
	)
	.sort(by = ['qseqid','sseqid','bitscore'], descending = [False,False,True])
	.unique(subset = ['qseqid', 'sseqid'], keep = 'first')
	.collect()

).write_csv('query_vs_neighboors.tsv', separator = '\t', include_header = True)
"