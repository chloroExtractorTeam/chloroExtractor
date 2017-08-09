chloroExtractor
===============

Table of contents
-----------------
1. Intoduction
2. Requirements
3. Installation
4. Usage
5. License


1. Introduction
---------------

The chloroExtractor is a perl based program which provides a pipeline for DNA extraction of chloroplast DNA from whole genom plant data.
Too huge amounts of chloroplast DNA can be a problem at the assambly of whole genom data, one solution for this problem can be a core extraction befor sequencing, but this can be expensive.
The chloroExtractor takes your whole genom data and extract the chloroplast DNA, so you can have your diffrent DNAs sepperated easily by the chloroExractor. Furthermore the chloroExtractor takes
the chloroplast DNA and trys to assamble it. This is possible because of the preserfed nature of the chloroplast primary but also the secondary structure.
Through k-mer filtering the k-mers which contains the chloroplast sequences get extracted and can then be used to assamble the chloroplast on a guided assamlby with several other chloroplasts. 



2. Requirements
---------------

*Jellyfish
*Spades

3. Installation
---------------

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

4. Usage
--------

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


5. License
----------

For License please refer to the LICENSE file