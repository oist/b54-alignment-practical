#!/bin/bash

grep "^UniRef90" result_drosophila_pax6_isoformA_blastp_uniref90.out | cut -d " " -f 1 | sort -u > blastp_uniref90_hits.list
grep "^NP" result_drosophila_pax6_isoformA_diamond_uniref90.out | cut -f 2 | sort -u > diamond_uniref90_hits.list
BLASTP_LINES=`wc -l blastp_uniref90_hits.list | cut -d " " -f 1`
DIAMOND_LINES=`wc -l diamond_uniref90_hits.list | cut -d " " -f 1`

COMMON_LINES=`grep -f diamond_uniref90_hits.list blastp_uniref90_hits.list | wc -l`

if [ $COMMON_LINES -eq $DIAMOND_LINES ]; then
	echo "All hits by diamond are present in blastp as well"
fi
