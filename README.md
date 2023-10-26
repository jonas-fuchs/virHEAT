[![DOI](https://zenodo.org/badge/639918477.svg)](https://zenodo.org/badge/latestdoi/639918477)
[![language](https://img.shields.io/badge/python-%3E3.9-green)](https://www.python.org/)
[![License: GPL v3](https://img.shields.io/github/license/jonas-fuchs/virheat)](https://www.gnu.org/licenses/gpl-3.0)
[![pypi version](https://img.shields.io/pypi/v/virheat)](https://pypi.org/project/virheat/)
[![pypi version](https://static.pepy.tech/badge/virheat)](https://pypi.org/project/virheat/)
[![CONDA](https://img.shields.io/conda/v/bioconda/virheat?label=conda%20version)](https://anaconda.org/bioconda/virheat)
[![CONDA](https://img.shields.io/conda/dn/bioconda/virheat?label=conda%20downloads)](https://anaconda.org/bioconda/virheat)

![Logo](./virheat.png)



**virHEAT is a tool to visualize vcfs as a heatmap and map mutations to respective genes.**



Ever wanted to have a condensed look at variant frequencies after mapping your raw reads to a viral/bacterial reference genome and compare multiple vcf files at the same time? Than virHEAT is for you. You can not only visualize the heatmap but also read in a gff3 file that lets you display genes harboring a mutation. This lightweight script was inspired by [snipit](https://github.com/aineniamh/snipit) and my [variant frequency plot](https://github.com/jonas-fuchs/SARS-CoV-2-analyses/tree/main/Heatmap), getting the best visualization features of both.

## SARS-CoV-2 example:

![Example](./example_data/example.png)

## Installation

### via pip (recommened):
```shell
pip install virheat
```
### via conda:
```shell
conda install -c bioconda virheat
```
### from this repo:
```shell
git clone https://github.com/jonas-fuchs/virHEAT
cd virHEAT
pip install -r requirements.txt
# or
pip install .
```
That was already it. To check if it worked:

```shell
virheat -v
```
You should see the current virHEAT version.

## Usage

```shell
usage: 	virheat <folder containing vcfs> <output dir> -l or -g [additional arguments]

```

**Arguments:**

```
positional arguments:
  input                 folder containing input files and output folder

options:
  -h, --help            show this help message and exit
  -l None, --genome-length None
                        length of the genome (needed if gff3 is not provided)
  -g None, --gff3-path None
                        path to gff3 (needed if length is not provided)
  -a [gene ...], --gff3-annotations [gene ...]
                        annotations to display from gff3 file (standard: gene). Multiple possible.
  -t 0, --threshold 0   display frequencies above this threshold (0-1)
  --delete, --no-delete
                        delete mutations that are present in all samples and their maximum frequency divergence is smaller than 0.5 (default: True)
  -n None, --delete-n None
                        do not show mutations that occur n times or less (default: Do not delete)                      
  --sort, --no-sort     sort sample names alphanumerically (default: False)
  --min-cov 20          display mutations covered at least x time (only if per base cov tsv files are provided)
  -v, --version         show program's version number and exit
```

You need to either provide the length of your reference genome or if you want to get the sequence annotation you will need to provide the gff3 file.

Additionally, you can also analyse if mutations are sufficiently covered and display non-covered cells in grey. For that first create a per base coverage tsv files for each bam file with [Qualimap](http://qualimap.conesalab.org/) and provide it in the same folder as the vcf files. Give them the same name as your vcf files.

<a href="https://www.buymeacoffee.com/jofox" target="_blank"><img src="https://www.buymeacoffee.com/assets/img/custom_images/orange_img.png" alt="Buy Me A Coffee" style="height: 41px !important;width: 174px !important;box-shadow: 0px 3px 2px 0px rgba(190, 190, 190, 0.5) !important;-webkit-box-shadow: 0px 3px 2px 0px rgba(190, 190, 190, 0.5) !important;" ></a>

---

**Important disclaimer:**
*The code is under the GPLv3 licence. The code is WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.*
