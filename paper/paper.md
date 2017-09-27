---
title: 'chloroExtractor: extraction and assembly of the chloroplast genome from whole genome shotgun data'
tags:
  - chloroplast
  - genome
  - assembly
authors:
 - name: Markus J Ankenbrand
   orcid: 0000-0002-6620-807X
   affiliation: 1,a
 - name: Simon Pfaff
   affiliation: 2,a
 - name: Niklas Terhoeven
   affiliation: 2,3
 - name: Musga Qureischi
   affiliation: 3,4
 - name: Maik Gündel
   orcid: 0000-0002-0502-4576
   affiliation: 3
 - name: Clemens L. Weiß
   affiliation: 5
 - name: Thomas Hackl
   orcid: 0000-0002-0022-320X
   affiliation: 6
 - name: Frank Förster
   orcid: 0000-0003-4166-5423
   affiliation: 2,3,7
affiliations:
 - name: Department of Animal Ecology and Tropical Biology (Zoology III), University of Würzburg, Germany
   index: 1
 - name: Center for Computational and Theoretical Biology, University of Würzburg
   index: 2
 - name: Department of Bioinformatics, University of Würzburg
   index: 3
 - name: Centre for Experimental Molecular Medicine, University Clinics Würzburg, Germany
   index: 4
 - name: Research Group for Ancient Genomics and Evolution, Department of Molecular Biology, Max Planck Institute for Developmental Biology, Tübingen, Germany
   index: 5
 - name: Department of Civil and Environmental Engineering, Massachusetts Institute of Technology
   index: 6
 - name: Fraunhofer Institute for Molecular Biology and Applied Ecology IME, Applied Ecology and Bioresources, Gießen, Germany
   index: 7
 - name: These authors contributed equally to this work
   index: a
date: 27 September 2017
bibliography: paper.bib
---

# Summary

This is an automated pipeline that extracts and reconstructs chloroplast genomes from whole genome shotgun data.
It is capable to assemble the incidental sequenced chloropast DNA, which is present in almost all plant sequencing projects, due to the extraction of whole cellular DNA.
It works by analyzing the k-mer distribution (determined with Jellyfish, [@marcais_fast_2011]) of the raw sequencing reads.
Usually the coverage of the chloroplast genome is much higher than that of the nuclear genome.
Using alignments to reference chloroplast sequences and the k-mer distribution candidate chloroplast reads are extracted from the complete set.
Afterwards, the targeted assembly of those sequences is much faster and yields less contigs compared to an assembly of all reads.
Assemblers usually fail to assemble chloroplast genomes as a single contig due to their structure, consisting of two single copy regions and an inverted repeat.
The size of the inverted repeat is in most cases multiple kilobasepairs in size, therefore it can not be resolved using short reads only.
However SPAdes [@Spades_2013] returns the assembly graph where the typical chloroplast structure can be recognized and reconstructed using the knowledge of its structure.
Using our demo set, one can achieve a single contig assembly of the chloroplast of *Spinacia oleracea* .
The final chloroplast sequence can be further annotated with tools like DOGMA [@dogma2004], cpGAVAS [@Liu2012] and VERDANT [@Mckain2017].
Such assemblies, can be used to remove chloroplast reads before a genomic assembly of the remaining nuclear DNA.
Moreover, chloroplast genomes are useful in phylogenetic reconstruction [@huang_2016] or barcoding applications [@Coissac_2016].
A similar tool, aiming the assembly of whole chloroplast genomes is the Python program [org.ASM](https://git.metabarcoding.org/org-asm/org-asm), but it is not production ready, yet.
Also plasmid SPAdes [@plasmidspades_2016] could possibly be used for this purpose although it is not intended for it.
In the future, we plan to use our chloroExtractor to screen NCBI's Sequence Read Archive [@sra2011] for chloroplast genomes in public sequencing datasets that are not yet available in chloroplast databases, eg. chloroDB [@chlordb2006] to broaden our knowledge about chloroplasts.

![Schematic workflow of chloroExtractor.](workflow.png)

# Acknowledgements

MJA was supported by a grant of the German Excellence Initiative to the Graduate School of Life Sciences, University of Würzburg.

# References

