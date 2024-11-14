"""
contains all plotting functions of virHEAT
"""

# BUILT-INS
import math
import numpy as np

# LIBS
import matplotlib.pyplot as plt
import matplotlib.patches as patches


def create_heatmap(ax, frequency_array, cmap):
    """
    create each cell of the heatmap
    """
    # cells
    y_start = 0
    for frequency_list in frequency_array:
        x_start = 0
        for freq in frequency_list:
            heatmap_cell = patches.Rectangle((x_start, y_start), 1, 1, alpha=1, edgecolor="grey", facecolor=cmap.to_rgba(freq))
            ax.add_patch(heatmap_cell)
            x_start += 1
        y_start += 1


def create_genome_vis(ax, genome_y_location, n_mutations, unique_mutations, start, stop):
    """
    create the genome rectangle, mutations and mappings to the heatmap
    """
    # colors
    mutation_type_colors = {
        "SNV": "dimgrey",
        "DEL": "red",
        "INS": "blue"
        }
    mutation_set = set()
    # define coordinates
    y_min = -genome_y_location
    y_max = -genome_y_location+genome_y_location/2

    # create blank rectangle for genome vis
    ax.add_patch(
        patches.FancyBboxPatch(
                            (0, y_min), n_mutations, y_max - y_min,
                            boxstyle="round,pad=-0.0040,rounding_size=0.03",
                            ec="none", fc="lightgrey"
                        )
    )
    # create mutation lines on the genome rectangle and the mapping to the respective cells
    x_start = 0
    length = stop - start
    for mutation in unique_mutations:
        mutation_attributes = mutation.split("_")
        mutation_color = mutation_type_colors[mutation_attributes[3]]
        mutation_set.add(mutation_attributes[3])
        mutation_x_location = n_mutations/length*(int(mutation_attributes[0])-start)
        # create mutation lines
        plt.vlines(x=mutation_x_location, ymin=y_min, ymax=y_max, color=mutation_color)
        # create polygon
        coordinates = [(x_start, 0), (x_start+1, 0), (mutation_x_location, y_max)]
        ax.add_patch(patches.Polygon(coordinates, facecolor=mutation_color, alpha=0.5))
        x_start += 1

    return mutation_set


def create_scores_vis(ax, genome_y_location, n_mutations, n_tracks, unique_mutations, start, stop, score_count, score_name):
    """
    create the scores rectangles, mappings to the reference
    """
    score_set = []
    y_zero = -genome_y_location - (n_tracks + score_count + 1) * (genome_y_location / 2)  # zero line
    length = stop - start

    # create list of tuples [(nt pos, score)]
    for mutation in unique_mutations:
        mutation_attributes = mutation.split("_")
        if not np.isnan(float(mutation_attributes[5])):
            score_set.append((int(mutation_attributes[0]), float(mutation_attributes[5])))
    # check if there is something to plot
    if score_set:
        # define normalization multiplier for the height of the score v lines
        max_value = max([abs(score[1]) for score in score_set])  # max abs in score set
        multiplier = max_value / abs((y_zero + genome_y_location / 4) - y_zero)
        # add text for scale and score_name
        if score_count % 2 == 0:
            score_name_loc, scale_loc, vline_loc = ax.get_xlim()[1] + 0.5, ax.get_xlim()[0] - 0.5, ax.get_xlim()[0]
            name_ha_aln, scale_ha_aln = "left", "right"
        else:
            score_name_loc, scale_loc, vline_loc = ax.get_xlim()[0] - 0.5, ax.get_xlim()[1] + 0.5, ax.get_xlim()[1]
            name_ha_aln, scale_ha_aln = "right", "left"
        # score name
        ax.text(score_name_loc, y_zero, score_name, ha=name_ha_aln, va='center')
        # scale
        ax.text(scale_loc, y_zero, "0", ha=scale_ha_aln, va='center')
        ax.text(scale_loc, y_zero - genome_y_location / 4, -max_value, ha=scale_ha_aln, va='center', color='red')
        ax.text(scale_loc, y_zero + genome_y_location / 4, max_value, ha=scale_ha_aln, va='center', color='blue')
        plt.vlines(x=vline_loc, ymin=y_zero - genome_y_location / 4, ymax=y_zero + genome_y_location / 4, color='black', linewidth=2)
        # add h line boundaries for scale
        if score_count is 1:
            plt.axhline(y=y_zero + genome_y_location / 4, color='grey', linestyle='--', linewidth=0.5)
        plt.axhline(y=y_zero - genome_y_location / 4, color='grey', linestyle='--', linewidth=0.5)
        # create zero line
        plt.axhline(y=y_zero, color='black', linestyle='-', linewidth=0.5)
        # create bars for scores
        for score in score_set:
            # define x value
            mutation_x_location = n_mutations / length * (score[0] - start)
            # create lines for score_set
            if score[1] < 0:
                plt.vlines(x=mutation_x_location, ymin=y_zero + score[1]/multiplier, ymax=y_zero, color="red", linestyle='-')
            else:
                plt.vlines(x=mutation_x_location, ymin=y_zero, ymax=y_zero + score[1] / multiplier, color="blue", linestyle='-')

        return True
    # if not the track is not created
    else:
        print("\033[31m\033[1mERROR:\033[0m Seems like there are no scores in the score set '{}' corresponding to the plotted mutation positions.".format(score_name))
        return False


def create_colorbar(threshold, cmap, min_y_location, n_samples, ax):
    """
    creates a custom colorbar and annotates the threshold
    """
    if n_samples >= 8:
        ticks = [0, 0.2, 0.4, 0.6, 0.8, 1]
        labels = [0, 0.2, 0.4, 0.6, 0.8, 1]
        if threshold + 0.1 in ticks or threshold - 0.1 in ticks:
            rounded_threshold = threshold
        else:
            rounded_threshold = round(threshold * 5) / 5
    else:
        ticks = [0, 0.5, 1]
        labels = [0, 0.5, 1]
        if threshold + 0.25 in ticks or threshold - 0.25 in ticks:
            rounded_threshold = threshold
        else:
            rounded_threshold = round(threshold * 2) / 2

    if rounded_threshold in ticks:
        ticks.remove(rounded_threshold)
        labels.remove(rounded_threshold)
    ticks.append(threshold)
    labels.append(f"threshold\n={threshold}")
    cbar = plt.colorbar(cmap, label="variant frequency", pad=0, shrink=n_samples/(min_y_location+n_samples), anchor=(0.1,1), aspect=15, ax=ax)
    cbar.set_ticks(ticks)
    cbar.set_ticklabels(labels)


def create_mutation_legend(mutation_set, min_y_location, n_samples, n_scores):
    """
    create a legend for the mutation type
    """
    legend_patches = []

    if "DEL" in mutation_set:
        legend_patches.append(patches.Patch(color="red", label='DEL'))
    if "INS" in mutation_set:
        legend_patches.append(patches.Patch(color="blue", label="INS"))
    if "SNV" in mutation_set:
        legend_patches.append(patches.Patch(color="dimgrey", label="SNV"))
    plt.legend(bbox_to_anchor=(1.01, 0.95-(n_samples/(min_y_location+n_samples+n_scores))), handles=legend_patches)


def create_axis(ax, n_mutations, min_y_location, n_samples, file_names, start, stop, genome_y_location, unique_mutations, reference_name):
    """
    create the axis of the plot
    """

    # define plot limits
    ax.set_xlim(0, n_mutations)
    ax.set_ylim(-min_y_location, n_samples)
    # define new ticks depending on the genome size
    axis_length = stop - start
    if n_mutations >= 20:
        xtick_dis = round(axis_length/6, -int(math.log10(axis_length / 6)) + 1)
        xtick_labels = [start, start + xtick_dis, start + xtick_dis*2, start + xtick_dis*3, start + xtick_dis*4, start + xtick_dis*5, stop]
    elif n_mutations >= 10:
        xtick_dis = round(axis_length / 3, -int(math.log10(axis_length / 3)) + 1)
        xtick_labels = [start, start + xtick_dis, start + xtick_dis * 2, stop]
    else:
        xtick_labels = [start, stop]
    xtick_labels = [int(tick) for tick in xtick_labels]
    # get the correct location of the genome pos on the axis
    xticks = [n_mutations/axis_length*(tick - start) for tick in xtick_labels]
    # set new ticks and change spines/yaxis
    ax.set_xticks(xticks, xtick_labels)
    # set y axis labels
    y_ticks = [idx+0.5 for idx in range(0, len(file_names))]
    y_ticks.append(-genome_y_location+genome_y_location/4)
    file_names.append(reference_name)
    ax.set_yticks(y_ticks, file_names)
    # set second x axis for mut pos
    secxtick_labels = []
    secxtick_ticks = []
    for idx, mut in enumerate(unique_mutations):
        mutation = mut.split("_")
        secxtick_labels.append(f"{mutation[1]}{mutation[0]}{mutation[2]}")
        secxtick_ticks.append(idx+0.5)
    secax = ax.secondary_xaxis("top")
    secax.set_xticks(secxtick_ticks, secxtick_labels, rotation=90)
    # set spine visibility
    ax.spines["bottom"].set_position(("data", -genome_y_location))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)


def create_gene_vis(ax, genes_with_mutations, n_mutations, y_size, n_tracks, start, stop, min_y_location, genome_y_location, colors_genes, n_scores, show_arrows = True):
    """
    create the vis for the gene
    """

    gene_annotations = []
    mult_factor = n_mutations/(stop-start)
    # define arrow head length for all arrows
    # dependent on how long the relative genes are to each other. If they are similar in length, it is roughly 1/8
    # of the shortest gene to display.
    if show_arrows:
        all_gene_lengths = [genes_with_mutations[x][0][1]-genes_with_mutations[x][0][0] for x in genes_with_mutations.keys()]
        if min(all_gene_lengths)*20 < max(all_gene_lengths):
            head_length = mult_factor*min(all_gene_lengths)
        elif min(all_gene_lengths)*10 < max(all_gene_lengths):
            head_length = mult_factor * min(all_gene_lengths) * 0.6
        elif min(all_gene_lengths)*5 < max(all_gene_lengths):
            head_length = mult_factor * min(all_gene_lengths) * 0.3
        else:
            head_length = mult_factor * min(all_gene_lengths) * 0.15

    for idx, gene in enumerate(genes_with_mutations):
        # define the plotting values for the patch
        if genes_with_mutations[gene][0][0] < start:
            x_value = 0
        else:
            x_value = mult_factor*(genes_with_mutations[gene][0][0]-start)
        y_value = -min_y_location+(n_tracks+n_scores-genes_with_mutations[gene][1])*genome_y_location/2
        if genes_with_mutations[gene][0][1] > stop:
            width = n_mutations - x_value
        else:
            width = mult_factor*(genes_with_mutations[gene][0][1]-start) - x_value
        height = genome_y_location/2
        # plot the patch
        if not show_arrows:
            ax.add_patch(
                patches.FancyBboxPatch(
                                    (x_value, y_value), width, height,
                                    boxstyle="round,pad=-0.0040,rounding_size=0.03",
                                    ec="black", fc=colors_genes[idx]
                                )
            )
        else:
            if genes_with_mutations[gene][0][2] == '-':
                x_value_arrow, arrow_width = x_value + width, -width
            else:
                x_value_arrow, arrow_width = x_value, width
            ax.add_patch(
                patches.FancyArrow(
                    x_value_arrow, y_value + height/2, dx=arrow_width,dy=0, width=height*0.7, head_width=height, head_length=head_length, length_includes_head=True,
                    ec="black", fc=colors_genes[idx]
                )
            )
        # define text pos for gene description inside or below the gene box, depending if it fits within
        if width > n_mutations/(y_size*8)*len(gene):
            gene_annotations.append(ax.text(x_value+width/2, y_value+height/2, gene, ha="center", va="center"))
        else:
            gene_annotations.append(ax.text(x_value+width/2, y_value-height/4, gene, rotation=40, rotation_mode="anchor", ha="right", va="bottom"))
