# chloroExtractor

## Introduction
The chloroExtractor is a perl based program which provides a pipeline for DNA extraction of chloroplast DNA from whole genome plant data.
Too huge amounts of chloroplast DNA can cast problems for the assembly of whole genome data.
One solution for this problem can be a core extraction before sequencing, but this can be expensive.
The chloroExtractor takes your whole genome data and extracts the chloroplast DNA, so you can have your different DNA separated easily by the chloroExractor.
Furthermore the chloroExtractor takes the chloroplast DNA and tries to assemble it.
This is possible because of the preserved nature of the chloroplasts primary and secondary structure.
Through k-mer filtering the k-mers which contain the chloroplast sequences get extracted and can then be used to assemble the chloroplast on a guided assembly with several other chloroplasts.



## Requirements
 - [Jellyfish](http://www.cbcb.umd.edu/software/jellyfish/ "Jellyfish K-mer counter")
 - [Spades](http://cab.spbu.ru/software/spades/ "SPAdes assamlber")
 - [Moose](http://search.cpan.org/~ether/Moose-2.2006/lib/Moose.pm "Moose Perl5-integration")
 - [Log4Perl](http://search.cpan.org/~mschilli/Log-Log4perl-1.49/lib/Log/Log4perl.pm "Log4Perl Perl5-Integration")
 - [Term::ProgressBar](http://search.cpan.org/~manwar/Term-ProgressBar-2.21/lib/Term/ProgressBar.pm "Term::ProgressBar Perl5-Integration")

## Installation
Just download the Archive and unpack it.

OR

Clone the directory form our github ->GITHUBLINK<-
```shell
git clone
```

make sure you uptade all the submodules  with:
```shell
cd chloroExtractor

git submodule init

git submodule clone --recursive
```

## Usage
To use the chloroExtractor, use the ptx executable in the bin/ folder


```shell

$ ./ptx --help

```

```shell
$ ./ptx [<OPTIONS>] -1 <FQ_1> -2 <FQ_2> -o <ID>

Options:
    -1|--reads
        Input reads file, first of pair.

    -2|--mates
        Input reads file, second of pair

    -d|--dir [ptx]
        Path to a working directory. Will be created. If exists, needs to be
        empty.

    -o|--out
        Output file basename, e.g. species id.

    --create-config
        Create a config file with default settings for user customization.

    -c|--config
        Use user customized config file. Superseeds default config.

    --continue=[TASKID TASKID ...] [TRUE]
        By default, the pipeline will check for a incomplete previous run
        and if possible continue after the last successful task of that run.
        Additionally you may provide task ids to specify a specific task -
        instead of the last task - to continue from.

    --redo [FALSE]
        Force pipeline to restart from the beginning, ignoring and
        overwriting previous results. Supersedes --continue.

    --stop-after=<TASKID>
        Stop the pipeline after the specified task.

    --skip=<TASKID/PATTERN TASKID/PATTERN ...>
        Skip specified tasks or tasks matching specified patterns (perl
        regex). If other tasks request results from skipped tasks, the
        pipeline will try to reuse results from previous runs. You need to
        take care, that these results still make sence in the current run.

    -V|--version
        Display version.

    -h|--help
        Display this help.
```

All the Options can and should be handled with the configuration file ptx.cfg, which is located in the mainfolder. With this config file you can handle the options for each step and task individual.
On default the chloroExtractor uses this config file, you can edit these one, or make your own one and uses the -c parameter to use it.

```shell

$ ./ptx -c ownptx.cfg -1 FQ_1 -2 FQ_2

```

## Data
The chloroExtractor uses unsortet Fastq files with paired end reads. Please make sure your reads are not sortet at all, otherwise there could be problems or even wrong results. 


## License
For License please refer to the LICENSE file
