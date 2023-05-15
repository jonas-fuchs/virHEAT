<img src="./virheat.png" alt="virHEAT logo" />


**virHEAT is a tool to visualize vcf files as a heatmap and map positions to a gff file.**



Ever wanted to have a condensed look at variant frequencies after mapping your raw reads to a viral/bacterial reference genome and compare multiple vcf files at the same time? Than virHEAT is for you.
You can not only visualize the heatmap but also read in a gff3 file that lets you display genes harboring a mutation.
This lightweight script was inspired by [snipit](https://github.com/aineniamh/snipit) and my [variant frequency plot](https://github.com/jonas-fuchs/SARS-CoV-2-analyses/tree/main/Heatmap), getting the best visualization features of both.

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

---

**Important disclaimer:**
*The code is under the GPLv3 licence. The code is WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.*
