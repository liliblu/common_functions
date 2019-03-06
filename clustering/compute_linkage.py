import matplotlib
matplotlib.use('agg')
import pandas as pd
import scipy.cluster.hierarchy
import scipy.spatial.distance.pdist
import argparse


if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Parse some arguments")
    parser.add_argument('--input_df', type=str)
    parser.add_argument('--output_prefix', type=str, default='clustered')

    args = parser.parse_args()

    input_df = pd.read_csv(args.input_df, sep='\t', index_col=0)
    output_prefix = args.output_prefix

    z = scipy.cluster.hierarchy.linkage(scipy.spatial.distance.pdist(input_df))

    z.to_csv('%s_linkage_matrix.txt' % output_prefix, sep='\t')
