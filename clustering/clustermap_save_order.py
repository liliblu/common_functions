import matplotlib
matplotlib.use('agg')
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

def fileToDict(tsv_map_file_name):
    map = dict()
    with open(tsv_map_file_name, 'r') as fh:
        for line in fh.readlines():
            map[line.split()[0]] = line.split()[1]
    return map


def makeClusterMapSaveOrder(df,
output_prefix,
sample_color_map,
group_color_map):

    group = df.columns
    column_colors = group.map(sample_color_map)

    group = df.index
    row_colors = group.map(sample_color_map)

    cg = sns.clustermap(df,
    cmap=cmap,
    vmin=-0.75,
    vmax=0.75,
    figsize=(40, 40),
    col_colors = column_colors,
    row_colors = row_colors,
    xticklabels=False,
    yticklabels=False,
    )

    cg.ax_row_dendrogram.set_visible(False)
    cg.ax_col_dendrogram.set_visible(False)

    row_order = cg.dendrogram_row.reordered_ind
    col_order = cg.dendrogram_col.reordered_ind

    with open('%s_row_order.txt' % output_prefix, 'w') as fh:
        for ind in row_order:
            fh.write('%s\n' %ind)

    with open('%s_col_order.txt' % output_prefix, 'w') as fh:
        for ind in col_order:
            fh.write('%s\n' %ind)

    handles = [mpatches.Patch(color=color, label=group) for group, color in group_color_map.iteritems()]

    fig = plt.gcf()
    fig.legend(handles=handles, bbox_to_anchor=(0.6, 0.0))

    cg.savefig('%s.png' % output_prefix, bbox_to_anchor='tight', pad_inches=1)

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Parse some arguments")
    parser.add_argument('--input_df', type=str)
    parser.add_argument('--output_prefix', type=str, default='clustered')
    parser.add_argument('--group_color_map', type=str, default='clustered')
    parser.add_argument('--sample_color_map', type=str, default='clustered')

    args = parser.parse_args()

    input_df = pd.read_csv(args.input_df, sep='\t', index_col=0)
    output_prefix = args.output_prefix

    sns.set(font = 'arial', style = 'white', color_codes=True, font_scale = 1.3)
    cmap = sns.diverging_palette(220, 10, as_cmap=True)
    cmap.set_bad('#BDBDBD')
    cmap.set_over('#9E031A')
    cmap.set_under('#0C4A60')

    group_color_map = fileToDict(args.group_color_map)
    sample_color_map = fileToDict(args.sample_color_map)

    makeClusterMapSaveOrder(input_df,
    output_prefix,
    sample_color_map,
    group_color_map)
