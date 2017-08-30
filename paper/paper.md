---
title: 'chloroExtractor: extraction and assembly of the chloroplast genome from whole genome shotgun data'
tags:
  - chloroplast
  - genome
  - assembly
authors:
 - name: Markus J Ankenbrand
   orcid: 0000-0002-6620-807X
   affiliation: 1
 - name: Simon Pfaff
   affiliation: 2
 - name: Niklas Terhoeven
   affiliation: 2
 - name: Musga
   affiliation: 3
 - name: Mike
   affiliation: 4
 - name: Clemens Weiß
   affiliation: 5
 - name: Thomas Hackl
   affiliation: 6
 - name: Frank Förster
   affiliation: 7
affiliations:
 - name: Department of Animal Ecology and Tropical Biology (Zoology III), University of Würzburg, Germany
   index: 1
date: 31 August 2017
bibliography: paper.bib
---

# Summary

This is an automated pipeline that extracts and reconstructs chloroplast genomes from whole genome shotgun data.
It works by analyzing the k-mer distribution (determined with jellyfish, cite) of the raw sequencing reads.
Usually the coverage of the chloroplast genome is much higher than that of the nuclear genome.
Using alignments to reference chloroplast sequences and the k-mer distribution candidate chloroplast reads are extracted and assembled.
This targeted assembly is much faster and yields less contigs compared to an assembly of all reads.
Assemblers usually fail to assemble chloroplast genomes as a single contig due to their structure, consisting of two single copy regions and an inverted repeat.
The size of the inverted repeat is in most cases multiple kilobasepairs in size, therefore it can not be resolved using short reads only.
However SPAdes (cite) returns the assembly graph where the typical chloroplast structure can be recognized and reconstructed using the knowledge of its structure.
This leads to a full chloroplast assembly as a single contig within minutes on this demo set: SRRxxxxxxx

Whole chloroplast genomes are useful in ...
The chloroExtractor bridges the gap between sequencing and annotation tools like DOGMA (cite), cpGAWAS (cite) and VERDANT (cite).
We plan to use the chloroExtractor also to screen SRA for chloroplast genomes in public datasets that are not yet available in chloroplast databases (cite NCBI and chloroDB).
The extraction of chloroplast sequences can be useful in order to remove chloroplast reads for the genomic assembly.

A similar tool to take note of is orgASM a dedicated organelle assembler written in python.
Also plasmid SPAdes could possibly be used for this purpose although it is not intended for it.

# References

