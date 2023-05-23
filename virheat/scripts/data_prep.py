"""
contains all data prep functions of virHEAT
"""

# BUILT-INS
import os
import re
import pathlib

# LIBS
import numpy as np


def get_vcf_files(path):
    """
    returns a list of vcf files in a paticular folder
    """
    return list(pathlib.Path(path).glob('*.vcf'))


def get_digit_and_alpha(filename):
    """
    get digits and alpha in file names
    """
    digit_match = re.search(r'\d+', filename)
    if digit_match:
        digit = digit_match.group()
        alpha = re.sub(r'\d+', '', filename)
    else:
        digit = ''
        alpha = filename
    return (alpha, int(digit) if digit else float('inf'))


def convert_string(string):
    """
    converts string to its right type
    """
    string = string.replace("\n", "")
    if string.isdecimal():
        return(int(string))
    elif string.replace('.', '', 1).isdecimal():
        return(float(string))
    else:
        return string


def read_vcf(vcf_file):
    """
    parse vcf files to dictionary
    """
    vcf_dict = {}

    # get header and values
    with open(vcf_file, "r") as f:
        header = [l.split("\t") for l in f if l.startswith('#CHROM')][0]
    # get each line as frequency_lists
    with open(vcf_file, "r") as f:
        lines = [l.split("\t") for l in f if not l.startswith('#')]
    # check if vcf is empty
    if not lines:
        print(f"WARNING: {vcf_file} is empty!")
    # get standard headers as keys
    for key in header[0:6]:
        vcf_dict[key] = []
    # functional effect
    vcf_dict["TYPE"] = []
    # info field
    for line in lines:
        for info in line[7].split(";"):
            if "=" in info:
                vcf_dict[info.split("=")[0]] = []
    # fill in dictionary
    for line in lines:
        # remember keys that have an entry already
        visited_keys = []
        for idx, key in enumerate(header[0:6]):
            vcf_dict[key].append(convert_string(line[idx]))
        # get mutation type
        if len(line[3]) == len(line[4]):
            vcf_dict["TYPE"].append("SNV")
        elif len(line[3]) < len(line[4]):
            vcf_dict["TYPE"].append("INS")
        elif len(line[3]) > len(line[4]):
            vcf_dict["TYPE"].append("DEL")
        visited_keys.extend(header[0:6])
        visited_keys.append("TYPE")
        # get data from info field
        for info in line[7].split(";"):
            if "=" in info:
                key, val = info.split("=")
                vcf_dict[key].append(convert_string(val))
                visited_keys.append(key)
        # append none for ech none vistited key
        for key in [k for k in vcf_dict.keys() if k not in visited_keys]:
            vcf_dict[key].append(None)

    return vcf_dict


def extract_vcf_data(vcf_files, threshold=0):
    """
    extract relevant vcf data
    """

    file_names = []
    frequency_lists = []

    for file in vcf_files:
        file_names.append(os.path.splitext(os.path.basename(file))[0])
        vcf_dict = read_vcf(file)
        frequency_list = []
        # write all mutation info in a '_' sep string
        for idx in range(0, len(vcf_dict["#CHROM"])):
            if not vcf_dict["AF"][idx] >= threshold:
                continue
            frequency_list.append(
                (f"{vcf_dict['POS'][idx]}_{vcf_dict['REF'][idx]}_{vcf_dict['ALT'][idx]}_{vcf_dict['TYPE'][idx]}", vcf_dict['AF'][idx])
            )
        frequency_lists.append(frequency_list)
    # sort by mutation index
    unique_mutations = sorted(
        {x[0] for li in frequency_lists for x in li}, key=lambda x: int(x.split("_")[0])
    )

    return vcf_dict["#CHROM"][0], frequency_lists, unique_mutations, file_names


def create_freq_array(unique_mutations, frequency_lists):
    """
    create an np array of the mutation frequencies
    """

    frequency_array = []

    for frequency_list in frequency_lists:
        frequencies = []
        for mutation in unique_mutations:
            af = [tup[1] for tup in frequency_list if tup[0] == mutation]
            if af:
                frequencies.append(af[0])
            else:
                frequencies.append(0)
        frequency_array.append(frequencies)

    return np.array(frequency_array)


def delete_common_mutations(frequency_array, unique_mutations):
    """
    delete rows of common mutations (non-zero) that are in the array
    """

    mut_to_del = []

    for idx in range(0, len(frequency_array[0])):
        for frequency_list in frequency_array:
            if frequency_list[idx] != 0:
                common_mut = True
            else:
                common_mut = False
                break
        if common_mut:
            mut_to_del.append(idx)

    for idx in sorted(mut_to_del, reverse=True):
        del unique_mutations[idx]

    return np.delete(frequency_array, mut_to_del, axis=1)


def parse_gff3(file):
    """
    parse gff3 to dictionary
    """

    gff3_dict = {}

    with open(file, "r") as gff3_file:
        for line in gff3_file:
            # ignore comments and last line
            if line.startswith("#") or line == "\n":
                continue
            gff_values = line.split("\t")
            gff3_ref_name = gff_values[0]
            # create keys
            if gff_values[2] not in gff3_dict:
                gff3_dict[gff_values[2]] = {}
            # parse the attribute line
            for attribute in gff_values[8].split(";"):
                identifier, val = attribute.split("=")
                # create a new dict for each ID
                if identifier == "ID" and identifier not in gff3_dict:
                    attribute_id = val
                    gff3_dict[gff_values[2]][attribute_id] = {}
                # add attributes
                if identifier != "ID":
                    gff3_dict[gff_values[2]][attribute_id][identifier] = val.replace("\n", "")
            # add start, stop and strand
            gff3_dict[gff_values[2]][attribute_id]["start"] = int(gff_values[3])
            gff3_dict[gff_values[2]][attribute_id]["stop"] = int(gff_values[4])
            gff3_dict[gff_values[2]][attribute_id]["strand"] = gff_values[6]

    gff3_file.close()

    return gff3_dict, gff3_ref_name


def get_genome_end(gff3_dict):
    """
    get the end of the genome from the region annotation
    """

    genome_end = 0

    for attribute in gff3_dict["region"].keys():
        stop = gff3_dict["region"][attribute]["stop"]
        if stop > genome_end:
            genome_end = stop

    return genome_end


def create_track_dict(unique_mutations, gff3_info):
    """
    create a dictionary of the genes that have mutations and assess in which
    track these genes should go in case they overlap
    """

    # find genes that have a mutation
    genes_with_mutations = set()

    for mutation in unique_mutations:
        # get the mutation from string
        mutation = int(mutation.split("_")[0])
        for gene_name in gff3_info["gene"]:
            if mutation in range(gff3_info["gene"][gene_name]["start"], gff3_info["gene"][gene_name]["stop"]):
                genes_with_mutations.add(
                    (gff3_info["gene"][gene_name]["gene"],
                     gff3_info["gene"][gene_name]["start"],
                     gff3_info["gene"][gene_name]["stop"],
                     gff3_info["gene"][gene_name]["strand"])
                )

    # create a dict and sort
    gene_dict = {element[0]: [element[1:4]] for element in genes_with_mutations}
    gene_dict = dict(sorted(gene_dict.items(), key=lambda x: x[1][0]))

    # remember for each track the largest stop
    track_stops = [0]

    for gene in gene_dict:
        track = 0
        # check if a start of a gene is smaller than the stop of the current track
        # -> move to new track
        while gene_dict[gene][0][0] < track_stops[track]:
            track += 1
            # if all prior tracks are potentially causing an overlap
            # create a new track and break
            if len(track_stops) <= track:
                track_stops.append(0)
                break
        # in the current track remember the stop of the current gene
        track_stops[track] = gene_dict[gene][0][1]
        # and indicate the track in the dict
        gene_dict[gene].append(track)

    return gene_dict, len(track_stops)
