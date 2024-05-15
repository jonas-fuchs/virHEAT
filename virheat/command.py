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
        "--name",
        type=str,
        metavar="virHEAT_plot.pdf",
        default="virHEAT_plot.pdf",
        help="plot name and file type (pdf, png, svg, jpg). Default: virHEAT_plot.pdf"
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
        "-z",
        "--zoom",
        type=int,
        metavar=("start", "stop"),
        nargs=2,
        help="restrict the plot to a specific genomic region."
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
        "-s", "--scores",
        metavar=('scores_file', 'pos_col', 'score_col', 'score_name'),
        nargs=4,
        action='append',
        help="specify scores to be added to the plot by providing a CSV file containing scores, along with its column for amino-acid positions, its column for scores, and descriptive score names (e.g., expression, binding, antibody escape, etc.). This option can be used multiple times to include multiple sets of scores."
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
    n_scores = 0
    if args.scores:
        reference_name, frequency_lists, unique_mutations, file_names = data_prep.extract_vcf_data(vcf_files, threshold=args.threshold, scores=1)
        n_scores = len(args.scores)
    else:
        reference_name, frequency_lists, unique_mutations, file_names = data_prep.extract_vcf_data(vcf_files, threshold=args.threshold)
    if args.zoom:
        unique_mutations = data_prep.zoom_to_genomic_regions(unique_mutations, args.zoom)
    frequency_array = data_prep.create_freq_array(unique_mutations, frequency_lists)

    # user specified delete options (removes mutations based on various rationales)
    if args.delete:
        frequency_array = data_prep.delete_common_mutations(frequency_array, unique_mutations)
    if args.delete_n is not None:
        frequency_array = data_prep.delete_n_mutations(frequency_array, unique_mutations, args.delete_n)

    # annotate low coverage if per base coverage from qualimap was provided
    data_prep.annotate_non_covered_regions(args.input[0], args.min_cov, frequency_array, file_names, unique_mutations)

    # define relative locations of all items in the plot
    n_samples, n_mutations = len(frequency_array), len(frequency_array[0])
    if n_mutations == 0:
        sys.exit("\033[31m\033[1mERROR:\033[0m Frequency array seems to be empty. There is nothing to plot.")
    if n_samples < 4:
        genome_y_location = 2
    else:
        genome_y_location = math.log2(n_samples)

    # gff3 data extraction and define y coordinates
    if args.gff3_path is not None and args.genome_length is not None:
        sys.exit("\033[31m\033[1mERROR:\033[0m Do not provide the -g and -l argument simultaneously!")
    elif args.gff3_path is not None:
        gff3_info, gff3_ref_name = data_prep.parse_gff3(args.gff3_path)
        # issue a warning if #CHROM and gff3 do not match
        if gff3_ref_name not in reference_name and reference_name not in gff3_ref_name:
            print("\033[31m\033[1mWARNING:\033[0m gff3 reference does not match the vcf reference!")
        genome_end = data_prep.get_genome_end(gff3_info)
        genes_with_mutations, n_tracks = data_prep.create_track_dict(unique_mutations, gff3_info, args.gff3_annotations)
    elif args.genome_length is not None:
        genome_end = args.genome_length
        n_tracks = 0
    else:
        sys.exit("\033[31m\033[1mERROR:\033[0m Provide either a gff3 file (-g) or the length (-l) of the genome which you used for mapping")

    # define min y coordinate
    if n_tracks != 0 or n_scores != 0:
        min_y_location = genome_y_location + genome_y_location / 2 * (n_tracks + n_scores + 1)
    else:
        min_y_location = genome_y_location

    # define size of the plot
    y_size = n_mutations*0.4
    if args.scores:
        x_size = y_size * (n_samples + min_y_location + len(args.scores)) / n_mutations
    else:
        x_size = y_size*(n_samples + min_y_location)/n_mutations
    x_size = x_size-x_size*0.15  # compensate of heatmap annotation

    # ini the fig
    fig, ax = plt.subplots(figsize=[y_size, x_size])

    # define boundaries for the plot
    if args.zoom:
        start, stop = args.zoom[0], args.zoom[1]
        # rescue plot if invalid zoom values are given
        if args.zoom[0] < 0:
            start = 0
        if args.zoom[1] > genome_end:
            stop = genome_end
    else:
        start, stop = 0, genome_end

    # plot all common elements
    cmap = cm.gist_heat_r
    cmap.set_bad('silver', 1.)
    cmap_cells = cm.ScalarMappable(norm=colors.Normalize(0, 1), cmap=cmap)
    plotting.create_heatmap(ax, frequency_array, cmap_cells)
    mutation_set = plotting.create_genome_vis(ax, genome_y_location, n_mutations, unique_mutations, start, stop)
    plotting.create_axis(ax, n_mutations, min_y_location, n_samples, file_names, start, stop, genome_y_location,
                         unique_mutations, reference_name, n_scores)
    plotting.create_mutation_legend(mutation_set, min_y_location, n_samples, n_scores)
    plotting.create_colorbar(args.threshold, cmap_cells, min_y_location, n_samples, ax, n_scores)
    # plot scores as track below the genome track
    if args.scores:
        score_count = 1
        for score_params in args.scores:
            score_count += 1
            scores_file, pos_col, score_col, score_name = score_params
            unique_scores = data_prep.extract_scores(unique_mutations, scores_file, pos_col, score_col)
            plotting.create_scores_vis(ax, genome_y_location, n_mutations, n_tracks, unique_scores, start, stop, score_name=score_name, score_count=score_count)



            #cmap_scores = plt.cm.get_cmap('coolwarm')  # blue to red colormap
            #if n_scoresets < 3:
            #    plotting.create_scores_cbar(cmap_scores, min_y_location, n_scoresets, score_set, score_name, score_count, ax)

        # plotting colorbars for scorsets if more than 3 scoresets - move colorbars under
        # if n_scoresets > 2:
        #     ax2 = fig.add_axes(ax.get_position(), frameon=False)
        #     ax2.set_position([ax.get_position().x0, ax.get_position().y0 - 1.1 * ax.get_position().height*0.5,
        #                       ax.get_position().width, 0.4])
        #     ax2.set_xticks([])
        #     ax2.set_yticks([])
        #     score_count = 0
        #     for score_params in args.scores:
        #         score_count += 1
        #         scores_file, pos_col, score_col, score_name = score_params
        #         unique_scores = data_prep.extract_scores(unique_mutations, scores_file, pos_col, score_col)
        #         score_set = plotting.create_scores_vis(ax, genome_y_location, n_mutations, unique_scores, start, stop, score_name=score_name, score_count=score_count, no_plot=1)
        #         cmap_scores = plt.cm.get_cmap('coolwarm')
        #         #plotting.create_scores_cbar(cmap_scores, min_y_location, n_scoresets, score_set, score_name, score_count, ax2)

    if args.gff3_path is not None:
        if genes_with_mutations:
            # distinct colors for the genes
            cmap_genes = plt.get_cmap('tab20', len(genes_with_mutations))
            colors_genes = [cmap_genes(i) for i in range(len(genes_with_mutations))]
            # plot gene track
            if args.scores and n_samples >= 4:
                min_y_location_sc = min_y_location - 1.7 * len(args.scores)
                plotting.create_gene_vis(ax, genes_with_mutations, n_mutations, y_size, n_tracks, start, stop, min_y_location_sc, genome_y_location, colors_genes)
            else:
                plotting.create_gene_vis(ax, genes_with_mutations, n_mutations, y_size, n_tracks, start, stop, min_y_location, genome_y_location, colors_genes)

    # create output folder
    if not os.path.exists(args.input[1]):
        os.makedirs(args.input[1])

    # save fig
    fig.savefig(os.path.join(args.input[1], args.name), bbox_inches="tight")

