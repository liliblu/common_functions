from __future__ import division
import pandas as pd, numpy as np, matplotlib.pyplot as plt, seaborn as sns
import matplotlib.patches as mpatches
import scipy.stats
import sys
import argparse
from datetime import datetime

sys.path.insert(0, "/Users/lili/Google Drive/Ruggles_lab/common_functions")
import commonFns

sns.set(font = 'arial', style = 'white', color_codes=True, font_scale = 1)
cmap = sns.cubehelix_palette(start=0.857, rot=0.00, gamma=1.5, hue=1, light=1, dark=0.2, reverse=False, as_cmap=True)
cmap.set_bad('#F5F5F5')

def makeHeatMap(heatmap_table, group_color_map, sample_color_map, output_prefix):
    group = heatmap_table.columns
    column_colors = group.map(sample_color_map)

    g = sns.clustermap(heatmap_table,
                           cmap=cmap,
                           col_cluster = False,
                           # row_cluster = False,
                           col_colors = column_colors,
                           xticklabels=False,
#                            yticklabels=False,
                        # standard_scale=0,
                        vmin=0,
#                         vmax=np.percentile(heatmap_table.values, 99.9),
                        vmax=30,
                           cbar_kws={'label':'# outliers'},
                          )
    g.ax_row_dendrogram.set_visible(False)
    plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    
    g.cax.set_position((0.15,0.12,0.03,0.6)) #move colorbar to right
    ax = g.ax_heatmap
    ax.set_ylabel('') #change the gene label
    
#     this chunk makes the legend the describes the different sample groups
    handles = [mpatches.Patch(color=color, label=group) for group, color in group_color_map.iteritems()]

    fig = plt.gcf()
    fig.legend(handles=handles, bbox_to_anchor=(0.6, 0.10))
    
    #save the plot
    plt.savefig('%s.pdf' %output_prefix, dpi=500, bbox_inches='tight', pad_inches=0.5)
    plt.show()
    plt.close()

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Input samples in two groups and df")
    parser.add_argument('--outliers_table', type=str,
                        help='tsv with proteins as rows and samples as columns. Phospho data also needs a "counts" column')
    parser.add_argument('--experiment_type', type=str,
                        choices=['phospho', 'not_phospho'],
                        default='phospho',
                        help='Phospho files need a "counts" column with total phosphosites measured per protein')
    parser.add_argument('--count_column_name', type=str,
                        default='counts',
                        help='Set for phospho data, if "counts" column is called something other than "counts"')
    parser.add_argument('--protein_column_name', type=str,
                        default='geneSymbol',
                        help='Specify column name for protein IDs')
    parser.add_argument('--fdr_cut_off', type=float,
                        default=0.05,
                        help='Set FDR cut off for which genes are considered sig diff')
    parser.add_argument('--output_prefix', type=str,
                        default=str(datetime.now().date()),
                        help='Default is current date. Set for prefix for gene lists and heatmap')
    parser.add_argument('--group1_label', type=str,
                        default='group1',
                        help='Label for group1 for heatmap')
    parser.add_argument('--group1_list', type=str,
                        help='List of samples from group 1. Samples separated by new lines, no header')
    parser.add_argument('--group2_label', type=str,
                        default='group2',
                        help='Label for group2 for heatmap')
    parser.add_argument('--group2_list', type=str,
                        help='List of samples from group 2. Samples separated by new lines, no header')
    parser.add_argument('--group_colors', type=str,
                        help='tsv, no header, with group names in 1st column and colors in 2nd column. Alternatively can map colors directly to samples. Group labels must match group labels given above.')

    args = parser.parse_args()

    outliers = pd.read_csv(args.outliers_table, sep='\t')
    experiment_type = args.experiment_type
    count_column_name = args.count_column_name
    
    if experiment_type == 'not_phospho':
        outliers[count_column_name] = 1
        
        
    protein_column_name = args.protein_column_name
    fdr_cut_off = args.fdr_cut_off

    output_prefix = args.output_prefix

    group1_label = args.group1_label
    group2_label = args.group2_label

    group1 = commonFns.fileToList(args.group1_list)
    group2 = commonFns.fileToList(args.group2_list)


# Assigning colors to samples
    if args.group_colors is not None:
        group_color_map = commonFns.fileToDict(group_colors)

        groups_dict = {sample:group1_label for sample in group1}
        groups_dict2 = {sample:group2_label for sample in group2}
        groups_dict.update(groups_dict2)

        sample_color_map = {sample:group_color_map[groups_dict[sample]] for sample in group1+group2}

    else:
        sample_color_map = {sample:'#571D41' for sample in group1}
        sample_color_map.update({sample:'#F5F5F5' for sample in group2})

        group_color_map = {group1_label:'#571D41', group2_label:'#F5F5F5'}


# Doing statistical test on different groups
    outliers['FDR'] = commonFns.testDifferentGroupsOutliers(group1,
                                                            group2,
                                                            outliers,
                                                            psite_count_column=count_column_name,
                                                            phospho=(experiment_type=='phospho'))

    outliers['significant'] = (outliers['FDR'] <= fdr_cut_off)
    sig_diff_count = sum(outliers['significant'])
    print('%s signficantly differential proteins' % sig_diff_count)

#If enough genes, make heatmap
    if sig_diff_count >= 2:
        heatmap_table = outliers.loc[(outliers['significant'] == True), [protein_column_name] + group1 + group2]
        heatmap_table = heatmap_table.set_index(protein_column_name)

        makeHeatMap(heatmap_table, group_color_map, sample_color_map, output_prefix)

#Write significantly different genes to a file
    if sig_diff_count > 0:
        up_in_group1 = outliers.loc[((outliers['significant']==True) & (outliers[group1].sum(axis=1) > outliers[group2].sum(axis=1))), protein_column_name]
        with open('%s.up_in_%s.txt' %(output_prefix, group1_label), 'w') as fh:
            for gene in up_in_group1:
                fh.write('%s\n'%gene)

        up_in_group2 = outliers.loc[((outliers['significant']==True) & (outliers[group1].sum(axis=1) < outliers[group2].sum(axis=1))), protein_column_name]
        with open('%s.up_in_%s.txt' %(output_prefix, group2_label), 'w') as fh:
            for gene in up_in_group2:
                fh.write('%s\n'%gene)
