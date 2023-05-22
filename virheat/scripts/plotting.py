"""
contains all plotting functions of virHEAT
"""

# BUILT-INS
import math

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


def create_genome_vis(ax, genome_y_location, n_mutations, unique_mutations, genome_end):
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
    # ax.add_patch(patches.Rectangle((0,y_min),n_mutations,y_max - y_min, alpha=1, edgecolor=None, facecolor="lightgrey"))
    ax.add_patch(
        patches.FancyBboxPatch(
                            (0, y_min), n_mutations, y_max - y_min,
                            boxstyle="round,pad=-0.0040,rounding_size=0.03",
                            ec="none", fc="lightgrey"
                        )
    )

    # create mutation lines on the genome rectangle and the mapping to the respective cells
    x_start = 0
    for mutation in unique_mutations:
        mutation_attributes = mutation.split("_")
        mutation_color = mutation_type_colors[mutation_attributes[3]]
        mutation_set.add(mutation_attributes[3])
        mutation_x_location = n_mutations/genome_end*int(mutation_attributes[0])
        # create mutation lines
        plt.vlines(x=mutation_x_location, ymin=y_min, ymax=y_max, color=mutation_color)
        # create polygon
        coordinates = [(x_start, 0), (x_start+1, 0), (mutation_x_location, y_max)]
        ax.add_patch(patches.Polygon(coordinates, facecolor=mutation_color, alpha=0.5))
        x_start += 1

    return mutation_set


def create_colorbar(threshold, cmap, min_y_location, n_samples):
    """
    creates a custom colorbar and annotates the threshold
    """
    ticks = [0, 0.2, 0.4, 0.6, 0.8, 1]
    labels = [0, 0.2, 0.4, 0.6, 0.8, 1]

    if threshold+0.1 in ticks or threshold-0.1 in ticks:
        rounded_threshold = threshold
    else:
        rounded_threshold = round(threshold*5)/5

    if rounded_threshold in ticks:
        ticks.remove(rounded_threshold)
        labels.remove(rounded_threshold)
    ticks.append(threshold)
    labels.append(f"threshold\n={threshold}")

    cbar = plt.colorbar(cmap, label="variant frequency", pad=0, shrink=n_samples/(min_y_location+n_samples), anchor=(0.1, 1))
    cbar.set_ticks(ticks)
    cbar.set_ticklabels(labels)


def create_mutation_legend(mutation_set, min_y_location, n_samples):
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

    plt.legend(bbox_to_anchor=(1.01, 0.95-(n_samples/(min_y_location+n_samples))), handles=legend_patches)


def create_axis(ax, n_mutations, min_y_location, n_samples, file_names, genome_end, genome_y_location, unique_mutations, reference_name):
    """
    create the axis of the plot
    """

    # define plot limits
    ax.set_xlim(0, n_mutations)
    ax.set_ylim(-min_y_location, n_samples)
    # define new ticks depending on the genome size
    xtick_dis = round(genome_end/6, -int(math.log10(genome_end/6))+1)
    xtick_labels = [0, xtick_dis, xtick_dis*2, xtick_dis*3, xtick_dis*4, xtick_dis*5, genome_end]
    xtick_labels = [int(tick) for tick in xtick_labels]
    # get the correct location of the genome pos on the axis
    xticks = [n_mutations/genome_end*tick for tick in xtick_labels]
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


def create_gene_vis(ax, genes_with_mutations, n_mutations, y_size, n_tracks, genome_end, min_y_location, genome_y_location, colors_genes):
    """
    create the vis for the gene
    """

    gene_annotations = []
    mult_factor = n_mutations/genome_end

    for idx, gene in enumerate(genes_with_mutations):
        start = (mult_factor*genes_with_mutations[gene][0][0], -min_y_location+(n_tracks-genes_with_mutations[gene][1])*genome_y_location/4-genome_y_location/4)
        stop = mult_factor*genes_with_mutations[gene][0][1]
        height = genome_y_location/2
        ax.add_patch(
            patches.FancyBboxPatch(
                                start, stop-start[0], height,
                                boxstyle="round,pad=-0.0040,rounding_size=0.03",
                                ec="black", fc=colors_genes[idx]
                            )
        )
        # define text pos for gene description inside or below the gene box, depending if it fits within
        if stop-start[0] > n_mutations/(y_size*8)*len(gene):
            gene_annotations.append(ax.text(start[0]+(stop-start[0])/2, start[1]+height/2, gene, ha="center", va="center"))
        else:
            gene_annotations.append(ax.text(start[0]+(stop-start[0])/2, start[1]-height/4, gene, ha="center", va="center"))
