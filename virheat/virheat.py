# BUILT-INS
import os
import pathlib

# LIBS
import numpy as np
import vcfpy as vcf
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import cm
from matplotlib import colors
from matplotlib import colormaps


def get_vcf_files(path):
    """
    returns a list of vcf files in a paticular folder
    """
    return list(pathlib.Path(path).glob('*.vcf'))


def extract_vcf_data(vcf_files, threshold=0):
    """
    extract relevant vcf data
    """

    file_names = []
    frequency_lists = []

    for file in vcf_files:
        file_names.append(os.path.splitext(os.path.basename(file))[0])
        vcf_reader = vcf.Reader.from_path(file)
        frequency_list = []
        # write all mutation info in a '_' sep string
        for line in vcf_reader:
            if not line.INFO["AF"] >= threshold:
                continue
            frequency_list.append(
                (f"{line.POS}_{line.REF}_{[alt.value for alt in line.ALT][0]}_{[alt.type for alt in line.ALT][0]}", line.INFO["AF"])
            )
        frequency_lists.append(frequency_list)
    # sort by mutation index
    unique_mutations = sorted(
        {x[0] for li in frequency_lists for x in li}, key=lambda x: int(x.split("_")[0])
    )

    return frequency_lists, unique_mutations, file_names


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

    return gff3_dict


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


def arg_parsing():
    """
    """
    print("todo")

def main():
    """
    """
    print("todo")


if __name__ == "__main__":
    main()
