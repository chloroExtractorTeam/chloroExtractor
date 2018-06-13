# chloroExtractor

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.883594.svg)](https://doi.org/10.5281/zenodo.883594)
[![Build Status](https://www.travis-ci.org/chloroExtractorTeam/chloroExtractor.svg)](https://www.travis-ci.org/chloroExtractorTeam/chloroExtractor)
[![Coverage Status](https://coveralls.io/repos/github/chloroExtractorTeam/chloroExtractor/badge.svg)](https://coveralls.io/github/chloroExtractorTeam/chloroExtractor)

## Introduction
The chloroExtractor is a perl based program which provides a pipeline for DNA extraction of chloroplast DNA from whole genome plant data.
Too huge amounts of chloroplast DNA can cast problems for the assembly of whole genome data.
One solution for this problem can be a core extraction before sequencing, but this can be expensive.
The chloroExtractor takes your whole genome data and extracts the chloroplast DNA, so you can have your different DNA separated easily by the chloroExractor.
Furthermore the chloroExtractor takes the chloroplast DNA and tries to assemble it.
This is possible because of the preserved nature of the chloroplasts primary and secondary structure.
Through k-mer filtering the k-mers which contain the chloroplast sequences get extracted and can then be used to assemble the chloroplast on a guided assembly with several other chloroplasts.



## Requirements
The version numbers given in parentheses are tested and known to work in this combination.
If you do a local install you can try to use other versions of some programs or modules but they are not guaranteed to work.
The docker container we provide will always contain a working combination of programs and modules.
### Required Software
 - [Jellyfish](http://www.cbcb.umd.edu/software/jellyfish/ "Jellyfish K-mer counter") (2.2.4)
 - [Spades](http://cab.spbu.ru/software/spades/ "SPAdes assamlber") (v3.10.1)
 - [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml "Bowtie2 Fast and sensitive read alignment") (2.2.6)
 - [NCBI-Blast+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download "BLAST (Basic Local Alignment Search Tool)") (2.2.31+)
 - [Samtools](http://www.htslib.org/ "Samtools Reading/writing/editing/indexing/viewing SAM/BAM/CRAM format") (0.1.19-96b5f2294a)
 - [Bedtools](http://bedtools.readthedocs.io/en/latest/ "bedtools: a powerful toolset for genome arithmetic") (v2.25.0)
 - [GNU R](https://www.r-project.org/ "The R Project for Statistical Computing") (3.2.3)
 - [Ghostscript](https://www.ghostscript.com/ "Ghostscript--an interpreter for the PostScript language and for PDF") (9.18)
 - [Python](https://www.python.org/ "www.python.org") (2.7.12)
 - [Perl](https://www.perl.org/ "www.perl.org") (v5.22.1)
### Required Perl modules
 - [Moose](http://search.cpan.org/~ether/Moose-2.2006/lib/Moose.pm "Moose Perl5-integration") (2.1604)
 - [Log::Log4Perl](http://search.cpan.org/~mschilli/Log-Log4perl-1.49/lib/Log/Log4perl.pm "Log4Perl Perl5-Integration") (1.44)
 - [Term::ProgressBar](http://search.cpan.org/~manwar/Term-ProgressBar-2.21/lib/Term/ProgressBar.pm "Term::ProgressBar Perl5-Integration") (2.17)
 - [Graph](http://search.cpan.org/dist/Graph/lib/Graph.pod "Graph - graph data structures and algorithms") (0.96)
 - [IPC::Run](http://search.cpan.org/~toddr/IPC-Run-0.96/lib/IPC/Run.pm "IPC::Run - system() and background procs w/ piping, redirs, ptys (Unix, Win32)") (0.94)
 - [File::Which](http://search.cpan.org/~plicease/File-Which-1.22/lib/File/Which.pm "File::Which - Perl implementation of the which utility as an API") (1.19)

## Installation
Install the requirements then clone the directory recursively
```shell
git clone --recursive https://github.com/chloroExtractorTeam/chloroExtractor
```

### Docker
[![Docker Automated build](https://img.shields.io/docker/automated/chloroextractorteam/chloroextractor.svg?style=plastic)](https://hub.docker.com/r/chloroextractorteam/chloroextractor/)
[![Docker Build Status](https://img.shields.io/docker/build/chloroextractorteam/chloroextractor.svg?style=plastic)](https://hub.docker.com/r/chloroextractorteam/chloroextractor/)

Our chloroExtractor is also available as a docker image.
Running chloroExtractor using that image requires the installation of docker and the permission to execute the `docker` commands.
Additionally, the docker container needs to be able to allocate enough memory (5GB are sufficient for the demo dataset).
In ubuntu memory for docker is usually not limited but on Mac OS X it is, refer to [this guide](https://stackoverflow.com/questions/32834082/how-to-increase-docker-machine-memory-mac/39720010#39720010) to increase the memory.
The data are mapped into the container as a volumne under `/data`.
Our chloroExtractor will run with `/data` as working directory.
Therefore, the output files will be stored inside the directory which was mapped into the container.
In case you are not using a user mapping, chloroExtractor will run with root priveleges and all created files will belong the `root` user.
For further information about docker and its security implications please visit their [website](https://docker.com).

```shell
docker pull chloroextractorteam/chloroextractor
docker run -v /location-of-input-data:/data --rm chloroextractorteam/chloroextractor -1 first_read.fq -2 second_read.fa [other options]
```

## Usage
To use the chloroExtractor, use the ptx executable in the bin/ folder

```shell
./ptx --help

```
or use the docker container:

``` script
docker run -v /location-of-input-data:/data --rm chloroextractorteam/chloroextractor --help
```
It returns a list of all mandatory parameters and optional setting.

```shell
$ ./ptx [<OPTIONS>] -1 <FQ_1> -2 <FQ_2> -d <OUTPUT-DIRECTORY>

Options:
    -1|--reads
        Input reads file, first of pair.

    -2|--mates
        Input reads file, second of pair

    -d|--dir [ptx]
        Path to a working directory. Will be created. If exists, needs to be
        empty.

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
## Input data
The chloroExtractor uses unsortet Fastq files with paired end reads. Please make sure your reads are not sortet at all, otherwise there could be problems or even wrong results.

## Example
An example data set can be downloaded from zenodo. As example we download the dataset into a folder and run chloroExtractor with the input files.

For preparation, a folder will be created and an example dataset will be downloaded:

```shell
# create a folder for the testrun, adjust this to your needs or use the current folder DATAFOLDER=$(pwd)
DATAFOLDER=/tmp/chloroExtractor-testrun
mkdir -p ${DATAFOLDER}
cd ${DATAFOLDER}

# download the example set and extract the sequencing reads
wget 'https://zenodo.org/record/884449/files/SRR5216995_1M.tar.bz2' -O - | tar xjf -
```
Afterwards, chloroExtractor can be run in command line mode:

```shell
# run chloroExtractor via command line (assuming all dependencies are installed and ptx folder is in PATH)
ptx -1 SRR5216995_1M_1.fastq -2 SRR5216995_1M_2.fastq
[17-09-21 13:42:42] [PipeWrap] Running ptx from the beginning, no previous runs detected.
[17-09-21 13:42:42] [PipeWrap] Running 'jf0': jellyfish count -t 8 -m 31 -s 500M -C -o jf0.jf /data/SRR5216995_1M_1.fastq /data/SRR5216995_1M_2.fastq
[...]
```

or using the docker container:

```shell
# other possibility is docker container based chloroExtractor (assuming that the user is allowed to run docker)
docker pull chloroextractorteam/chloroextractor # ensure the latest version from docker hub
# this binds the DATAFOLDER from above into the docker container you can also use the path directly instead of the variable
docker run -v ${DATAFOLDER}:/data --rm chloroextractorteam/chloroextractor -1 SRR5216995_1M_1.fastq -2 SRR5216995_1M_2.fastq
[17-09-21 13:52:30] [PipeWrap] Running ptx from the beginning, no previous runs detected.
[17-09-21 13:52:30] [PipeWrap] Running 'jf0': jellyfish count -t 8 -m 31 -s 500M -C -o jf0.jf /data/SRR5216995_1M_1.fastq /data/SRR5216995_1M_2.fastq
[...]
```

Both runs result in a final chloroplast assembly in the file `fcg.fa`.

Another more detailed example is available at our [demo](DEMO.md).

## Changelog
Version 1.0.4 is archived as [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1288679.svg)](https://doi.org/10.5281/zenodo.1288679). It updates the fastg-parser to version v0.6.0 and therefore fixes #101

Version 1.0.3 is archived as [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1179297.svg)](https://doi.org/10.5281/zenodo.1179297). It includes a test set and a patch for the divede by zero bug.

Version 1.0.2 is archived as [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1148955.svg)](https://doi.org/10.5281/zenodo.1148955) and was created after review process in [The Journal of Open Source Software](http://joss.theoj.org/), we added Daniel Amsel to the acknowledgement section.

Version 1.0.1 is archived as [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1147434.svg)](https://doi.org/10.5281/zenodo.1147434) and was created after review process in [The Journal of Open Source Software](http://joss.theoj.org/)

Version 1.0.0 is archived as [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.998262.svg)](https://doi.org/10.5281/zenodo.998262) and used for submission to [The Journal of Open Source Software](http://joss.theoj.org/)

## License
For License please refer to the LICENSE file
