#!/bin/bash

INPUTFILE="contigs_renamed.fa.contigs"
READS="../../cutoff_reads_out.fq"
MATES="../../cutoff_mates_out.fq"

OUT_READS="cutoff_reads_mapped.fq"
OUT_MATES="cutoff_mates_mapped.fq"

##
## mapping parameter
##

MAPPINGREGION=200            # how long are the regions onto which will be mapped
MINCONTIGLENGTH=500          # minimal contig length

##
## filenames used for temporary storage
##

COORDSFILE="$INPUTFILE"_end.coords
CONTIGENDS="$INPUTFILE"_end.fa

BOWTIE2_BASE="$INPUTFILE"_end
BOWTIE2_ALIGNED="$BOWTIE2_BASE"_al
BOWTIE2_UNALIGNED="$BOWTIE2_BASE"_un
BOWTIE2_SAM="$BOWTIE2_BASE".sam

MAPPED_READ_IDS="$INPUTFILE".al.ids

grep ">"  "$INPUTFILE" | perl -ne 's/>//;chomp; print "$_\t0\t$MAPPINGREGION\n$_\t-$MAPPINGREGION\t$MAPPINGREGION\n"' > "$COORDSFILE"
SeqFilter --in "$INPUTFILE" --min-length "$MINCONTIGLENGTH" --out - | SeqFilter --substr "$COORDSFILE" --substr-perl-style --out $CONTIG_ENDS
/storage/software/amd64/bowtie2/bowtie2-build "$CONTIG_ENDS" "$BOWTIE2BASE"
/storage/software/amd64/bowtie2/bowtie2 --phred33 --al-conc "$BOWTIE2_ALIGNED" --un-conc "BOWTIE2_UNALIGNED" --no-unal -x "$BOWTIE2_BASE" -1 "$READS" -2 "$MATES" -S "$BOWTIE2_SAM"
grep -v "^@" "$BOWTIE2_SAM" | cut -f1 >"$MAPPED_READ_IDS"
SeqFilter --in "$READS" --ids "$MAPPED_READ_IDS" --out "$OUT_READS"
SeqFilter --in "$MATES" --ids "$MAPPED_READ_IDS" --out "$OUT_MATES"


