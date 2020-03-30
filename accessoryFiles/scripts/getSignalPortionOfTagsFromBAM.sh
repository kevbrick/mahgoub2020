#!/bin/bash

bam=$1
bed=$2
idx=$3

if test -f $idx; then 
  >&2 echo "IDX: $idx"
else
  idx='mm10_genome.fa.fai'
  >&2 echo "IDX: $idx"
fi

sbed=${bed/bed/gS.bed}

>&2 echo STDERR "sortBed -i $bed -g $idx |cut -f1-3 >$sbed"
>&2 echo STDERR $idx
>&2 echo STDERR "intersectBed -a $bam -b $sbed -c -bed -sorted -g $idx |cut -f13"
>&2 echo STDERR "####### START ########"

tmp=$RAND".bed"
grep -P '^chr[0123456789XY]+\s' $bed |cut -f1-3 >$tmp
sortBed -i $tmp -g $idx >$sbed 

intersectBed -a $bam -b $sbed -c -bed -sorted -g $idx |cut -f13| perl -lane 'if($F[0] > 0){$in++}else{$out++}; $tot=$in+$out; print join("\t","'$bam'","'$sbed'",$in,$out,$in/$tot*100,$out/$tot*100)' |tail -n1

