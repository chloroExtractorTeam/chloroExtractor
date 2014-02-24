#!/bin/bash

grep ">" contigs_renamed.fa.contigs | perl -ne 's/>//;chomp; print "$_\t0\t200\n$_\t-200\t200\n"' >contigs_renamed.fa.contigs_end.coords
SeqFilter --in contigs_renamed.fa.contigs --min-length 500 --out - | SeqFilter --substr contigs_renamed.fa.contigs_end.coords --substr-perl-style --out contigs_renamed.fa.contigs_end.fa
/storage/software/amd64/bowtie2/bowtie2-build contigs_renamed.fa.contigs_end.fa contigs_renamed.fa.contigs_end
/storage/software/amd64/bowtie2/bowtie2 --phred33 --al-conc contigs_renamed.fa.contigs_end_al --un-conc contigs_renamed.fa.contigs_end_un --no-unal -x contigs_renamed.fa.contigs_end -1 ../../cutoff_reads_out.fq -2 ../../cutoff_mates_out.fq -S contigs_renamed.fa.contigs_end.sam 
grep -v "^@" contigs_renamed.fa.contigs_end.sam | cut -f1 >contigs_renamed.fa.contigs_end.al.ids
SeqFilter --in ../../cutoff_reads_out.fq --ids contigs_renamed.fa.contigs_end.al.ids --out cutoff_reads_mapped.fq
SeqFilter --in ../../cutoff_mates_out.fq --ids contigs_renamed.fa.contigs_end.al.ids --out cutoff_mates_mapped.fq

