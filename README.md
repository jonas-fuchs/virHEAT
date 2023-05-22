[![CodeFactor](https://www.codefactor.io/repository/github/jonas-fuchs/virheat/badge)](https://www.codefactor.io/repository/github/jonas-fuchs/virheat)
[![DOI](https://zenodo.org/badge/639918477.svg)](https://zenodo.org/badge/latestdoi/639918477)
[![language](https://img.shields.io/badge/python-%3E3.9-green)](https://www.python.org/)
[![License: GPL v3](https://img.shields.io/github/license/jonas-fuchs/virheat)](https://www.gnu.org/licenses/gpl-3.0)
[![pypi version](https://img.shields.io/pypi/v/virheat)](https://pypi.org/project/virheat/)

[![Logo](https://github.com/jonas-fuchs/virHEAT/blob/master/virheat.png)](https://github.com/jonas-fuchs/virHEAT/blob/master/virheat.png)



**virHEAT is a tool to visualize vcfs as a heatmap and map mutations to respective genes.**



Ever wanted to have a condensed look at variant frequencies after mapping your raw reads to a viral/bacterial reference genome and compare multiple vcf files at the same time? Than virHEAT is for you. You can not only visualize the heatmap but also read in a gff3 file that lets you display genes harboring a mutation. This lightweight script was inspired by [snipit](https://github.com/aineniamh/snipit) and my [variant frequency plot](https://github.com/jonas-fuchs/SARS-CoV-2-analyses/tree/main/Heatmap), getting the best visualization features of both.

## SARS-CoV-2 example:

[![Example](https://github.com/jonas-fuchs/virHEAT/blob/master/example.png)](https://github.com/jonas-fuchs/virHEAT/blob/master/example.png)

## Installation

### via pip (recommened):
```shell
pip install virheat
```
### from this repo:
```shell
git clone https://github.com/jonas-fuchs/varVAMP
cd virHEAT
```
and then install virHEAT with:
```shell
pip install -r requirements.txt
```
or:
```shell
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
  input                 folder containing vcf files and output folder

optional arguments:
  -h, --help            show this help message and exit
  -l None, --genome-length None
                        length of the genome (needed if gff3 is not provided)
  -g None, --gff3-path None
                        path to gff3 (needed if length is not provided)
  -t 0, --threshold 0   display frequencies above this threshold
  --delete, --no-delete
                        delete mutations with frequencies present in all
                        samples (default: True)
  --sort, --no-sort     sort alphanumerically (default: False)
  -v, --version         show program's version number and exit
```

You need to either provide the length of your reference genome or if you want to get the sequence annotation you will need to provide the gff3 file.

## Planned features

- Clustering
- Interactive mode

---

**Important disclaimer:**
*The code is under the GPLv3 licence. The code is WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.*
