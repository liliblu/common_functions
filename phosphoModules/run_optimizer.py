import argparse
import sys
sys.path.insert(0, '/gpfs/data/ruggleslab/home/lmb529/phos_corr/methods/phosphoModules')
import clusterOptimize

gctFile = '/gpfs/home/lmb529/ruggleslabHome/phos_corr/brca-prospective/normed/prosp-brca-v3.0-phosphoproteome-unfiltered-dedup-prot-normed.txt'
samplesFile = '/gpfs/home/lmb529/ruggleslabHome/phos_corr/brca-prospective/normed/sample_roster.txt'
gmtFile = '/gpfs/home/lmb529/ruggleslabHome/phos_corr/positive_controls/ptm.sig.db.all.flanking.human.v1.9.0.gmt'
writeToFilePrefix = '/gpfs/home/lmb529/ruggleslabHome/phos_corr/brca-prospective/normed/brca_pros_normed'

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Parse some arguments")
    parser.add_argument('--std', nargs='+', type=float)
    parser.add_argument('--nmembers', nargs='+', type=int)

    args = parser.parse_args()

    std = args.std
    nmem = args.nmembers

    scores = clusterOptimize.optimize(gctFile,
                              samplesFile,
                              gmt_file=gmtFile,
                              writeToFilePrefix=writeToFilePrefix,
                              gene_col='accession_number',
                              psite_col='variableSites',
                              std_range=std,
                              nmembers_range=nmem,
                              write_labels=True)

    output_file = '/gpfs/home/lmb529/ruggleslabHome/phos_corr/brca-prospective/optimization_%s_std.txt' %('.'.join([str(x) for x in std]))
    scores.to_csv(output_file, index=False, sep='\t')
