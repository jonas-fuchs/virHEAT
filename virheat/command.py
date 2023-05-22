"""
contains the main workflow of virHEAT
"""

# BUILT-INS
import os
import sys
import argparse

# LIB
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import colors
from matplotlib import colormaps

# virHEAT
from virheat.scripts import data_prep
from virheat.scripts import plotting
from virheat import __version__
from . import _program


def get_args(sysargs):
    """
    arg parsing for virheat
    """
    parser = argparse.ArgumentParser(
        prog=_program,
        usage='''\tvirheat <folder containing vcfs> <output dir> -l or -g [additional arguments]''')
    parser.add_argument(
        "input",
        nargs=2,
        help="folder containing vcf files and output folder"
    )
    parser.add_argument(
        "-l",
        "--genome-length",
        type=int,
        metavar="None",
        default=None,
        help="length of the genome (needed if gff3 is not provided)"
    )
    parser.add_argument(
        "-g",
        "--gff3-path",
        type=str,
        metavar="None",
        default=None,
        help="path to gff3 (needed if length is not provided)"
    )
    parser.add_argument(
        "-t",
        "--threshold",
        type=float,
        metavar="0",
        default=0,
        help="display frequencies above this threshold"
    )
    parser.add_argument(
        "--delete",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="delete mutations with frequencies present in all samples"
    )
    parser.add_argument(
        "--sort",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="sort alphanumerically"
    )
    parser.add_argument(
        "-v",
        "--version",
        action='version',
        version=f"virheat {__version__}"
    )
    if len(sysargs) < 1:
        parser.print_help()
        sys.exit(-1)
    else:
        return parser.parse_args(sysargs)


def main(sysargs=sys.argv[1:]):
    """
    main function for data extraction and plotting
    """

    # parse args
    args = get_args(sysargs)

    # get vcf files and sort
    vcf_files = data_prep.get_vcf_files(args.input[0])
    if args.sort:
        vcf_files = sorted(vcf_files, key=lambda x: data_prep.get_digit_and_alpha(os.path.basename(x)))

    # extract vcf info
    reference_name, frequency_lists, unique_mutations, file_names = data_prep.extract_vcf_data(vcf_files, threshold=args.threshold)
    frequency_array = data_prep.create_freq_array(unique_mutations, frequency_lists)
    if args.delete:
        frequency_array = data_prep.delete_common_mutations(frequency_array, unique_mutations)

    # define relative locations of all items in the plot
    n_samples = len(frequency_array)
    n_mutations = len(frequency_array[0])
    genome_y_location = n_samples/4

    # gff3 data extraction
    if args.gff3_path is not None and args.genome_length is not None:
        sys.exit("\033[31m\033[1mERROR:\033[0m Do not provide the -g and -l argument simultaneously!")
    elif args.gff3_path is not None:
        gff3_info, gff3_ref_name = data_prep.parse_gff3(args.gff3_path)
        # issue a warning if #CHROM and gff3 do not match
        if gff3_ref_name not in reference_name and reference_name not in gff3_ref_name:
            print("\033[31m\033[WARNING:\033[0m gff3 reference does not match the vcf reference!")
        genome_end = data_prep.get_genome_end(gff3_info)
        genes_with_mutations, n_tracks = data_prep.create_track_dict(unique_mutations, gff3_info)
        # define space for the genome vis tracks
        min_y_location = genome_y_location + genome_y_location/4 * (n_tracks+2)
    elif args.genome_length is not None:
        genome_end = args.genome_length
        min_y_location = genome_y_location
    else:
        sys.exit("\033[31m\033[1mERROR:\033[0m Provide either a gff3 file (-g) or the length (-l) of the genome which you used for mapping")

    # define size of the plot
    y_size = (n_mutations)*0.4
    x_size = y_size*(n_samples+min_y_location)/n_mutations
    x_size = x_size-x_size*1/6  # compensate of heatmap annotation

    # ini the fig
    fig, ax = plt.subplots(figsize=[y_size, x_size])

    # plot all elements
    cmap_cells = cm.ScalarMappable(norm=colors.Normalize(0, 1), cmap=colormaps["gist_heat_r"])
    plotting.create_heatmap(ax, frequency_array, cmap_cells)
    mutation_set = plotting.create_genome_vis(ax, genome_y_location, n_mutations, unique_mutations, genome_end)
    if args.gff3_path is not None:
        # distinct colors for the genes
        cmap_genes = plt.get_cmap("tab20")
        colors_genes = [cmap_genes(i) for i in range(len(genes_with_mutations))]
        # plot gene track
        plotting.create_gene_vis(ax, genes_with_mutations, n_mutations, y_size, n_tracks, genome_end, min_y_location, genome_y_location, colors_genes)
    plotting.create_axis(ax, n_mutations, min_y_location, n_samples, file_names, genome_end, genome_y_location, unique_mutations, reference_name)
    plotting.create_colorbar(args.threshold, cmap_cells, min_y_location, n_samples)
    plotting.create_mutation_legend(mutation_set, min_y_location, n_samples)

    # create output folder
    if not os.path.exists(args.input[1]):
        os.makedirs(args.input[1])

    # save fig
    fig.savefig(os.path.join(args.input[1], "virHEAT_plot.pdf"), bbox_inches="tight")
