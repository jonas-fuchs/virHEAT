"""
contains the main workflow of virHEAT
"""
# BUILT-INS
import os
import sys
import math
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
        usage='''\tvirheat <folder containing input files (vcf/tsv)> <output dir> -l or -g [additional arguments]''')
    parser.add_argument(
        "input",
        nargs=2,
        help="folder containing input files and output folder"
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
        "-a",
        "--gff3-annotations",
        type=str,
        action="store",
        metavar="gene",
        nargs="*",
        default=["gene"],
        help="annotations to display from gff3 file (standard: gene). Multiple possible."
    )
    parser.add_argument(
        "-t",
        "--threshold",
        type=float,
        metavar="0",
        default=0,
        help="display frequencies above this threshold (0-1)"
    )
    parser.add_argument(
        "--delete",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="delete mutations that are present in all samples and their maximum frequency divergence is smaller than 0.5"
    )
    parser.add_argument(
        "-n",
        "--delete-n",
        type=int,
        metavar='None',
        default=None,
        help="do not show mutations that occur n times or less (default: Do not delete)"
    )
    parser.add_argument(
        "--sort",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="sort sample names alphanumerically"
    )
    parser.add_argument(
        "--min-cov",
        type=int,
        metavar="20",
        default=20,
        help="display mutations covered at least x time (only if per base cov tsv files are provided)"
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
    vcf_files = data_prep.get_files(args.input[0], "vcf")
    if args.sort:
        vcf_files = sorted(vcf_files, key=lambda x: data_prep.get_digit_and_alpha(os.path.basename(x)))

    # extract vcf info
    reference_name, frequency_lists, unique_mutations, file_names = data_prep.extract_vcf_data(vcf_files, threshold=args.threshold)
    frequency_array = data_prep.create_freq_array(unique_mutations, frequency_lists)
    # user specified delete options (removes mutations based on various rationales)
    if args.delete:
        frequency_array = data_prep.delete_common_mutations(frequency_array, unique_mutations)
    if args.delete_n is not None:
        frequency_array = data_prep.delete_n_mutations(frequency_array, unique_mutations, args.delete_n)
    # annotate low coverage if per base coveage from qualimap was provided
    data_prep.annotate_non_covered_regions(args.input[0], args.min_cov, frequency_array, file_names, unique_mutations)

    # define relative locations of all items in the plot
    n_samples = len(frequency_array)
    n_mutations = len(frequency_array[0])
    if n_mutations == 0:
        sys.exit("\033[31m\033[1mERROR:\033[0m Frequency array seems to be empty. There is nothing to plot.")
    if n_samples < 4:
        genome_y_location = 2
    else:
        genome_y_location = math.log2(n_samples)

    # gff3 data extraction
    if args.gff3_path is not None and args.genome_length is not None:
        sys.exit("\033[31m\033[1mERROR:\033[0m Do not provide the -g and -l argument simultaneously!")
    elif args.gff3_path is not None:
        gff3_info, gff3_ref_name = data_prep.parse_gff3(args.gff3_path)
        # issue a warning if #CHROM and gff3 do not match
        if gff3_ref_name not in reference_name and reference_name not in gff3_ref_name:
            print("\033[31m\033[1mWARNING:\033[0m gff3 reference does not match the vcf reference!")
        genome_end = data_prep.get_genome_end(gff3_info)
        genes_with_mutations, n_tracks = data_prep.create_track_dict(unique_mutations, gff3_info, args.gff3_annotations)
        # define space for the genome vis tracks
        min_y_location = genome_y_location + genome_y_location/2 * (n_tracks+1)
    elif args.genome_length is not None:
        genome_end = args.genome_length
        min_y_location = genome_y_location
    else:
        sys.exit("\033[31m\033[1mERROR:\033[0m Provide either a gff3 file (-g) or the length (-l) of the genome which you used for mapping")

    # define size of the plot
    y_size = n_mutations*0.4
    x_size = y_size*(n_samples+min_y_location)/n_mutations
    x_size = x_size-x_size*0.15  # compensate of heatmap annotation

    # ini the fig
    fig, ax = plt.subplots(figsize=[y_size, x_size])

    # plot all elements
    cmap = cm.gist_heat_r
    cmap.set_bad('silver', 1.)
    cmap_cells = cm.ScalarMappable(norm=colors.Normalize(0, 1), cmap=cmap)
    plotting.create_heatmap(ax, frequency_array, cmap_cells)
    mutation_set = plotting.create_genome_vis(ax, genome_y_location, n_mutations, unique_mutations, genome_end)
    if args.gff3_path is not None:
        if genes_with_mutations:
            # distinct colors for the genes
            cmap_genes = plt.get_cmap('tab20', len(genes_with_mutations))
            colors_genes = [cmap_genes(i) for i in range(len(genes_with_mutations))]
            # plot gene track
            plotting.create_gene_vis(ax, genes_with_mutations, n_mutations, y_size, n_tracks, genome_end, min_y_location, genome_y_location, colors_genes)
    plotting.create_axis(ax, n_mutations, min_y_location, n_samples, file_names, genome_end, genome_y_location, unique_mutations, reference_name)
    plotting.create_colorbar(args.threshold, cmap_cells, min_y_location, n_samples, ax)
    plotting.create_mutation_legend(mutation_set, min_y_location, n_samples)

    # create output folder
    if not os.path.exists(args.input[1]):
        os.makedirs(args.input[1])

    # save fig
    fig.savefig(os.path.join(args.input[1], "virHEAT_plot.pdf"), bbox_inches="tight")

