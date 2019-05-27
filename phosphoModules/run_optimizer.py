import argparse
import sys
sys.path.insert(0, '/gpfs/data/ruggleslab/home/lmb529/phos_corr/methods/phosphoModules')
import testAll

gctFile = '/gpfs/home/lmb529/ruggleslabHome/phos_corr/brca-prospective/brca_prosp_v2.1_phosphoproteome-ratio-norm-NArm.gct_site-centric-seqwin_localized_n130x27352.gct'
samplesFile = '/gpfs/home/lmb529/ruggleslabHome/phos_corr/brca-prospective/sample_roster.txt'
gmtFile = '/gpfs/home/lmb529/ruggleslabHome/phos_corr/positive_controls/ptm.sig.db.all.flanking.human.v1.9.0.gmt'
writeToFilePrefix = '/gpfs/home/lmb529/ruggleslabHome/phos_corr/brca-prospective/brca_pros'

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Parse some arguments")
    parser.add_argument('--std', nargs='+', type=float)
    parser.add_argument('--nmembers', nargs='+', type=int)

    args = parser.parse_args()

    std = args.std
    nmem = args.nmembers

    scores = testAll.optimize(gctFile,
                              samplesFile,
                              gmt_file=gmtFile,
                              writeToFilePrefix=writeToFilePrefix,
                              std_range=std,
                              nmembers_range=nmem,
                              write_labels=True)

    output_file = '/gpfs/home/lmb529/ruggleslabHome/phos_corr/brca-prospective/optimization_%s_std.txt' %('.'.join([str(x) for x in std]))
    scores.to_csv(output_file, index=False, sep='\t')
