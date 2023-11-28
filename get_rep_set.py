import qiime2
import biom
import pandas as pd
import sys

t = qiime2.Artifact.load(sys.argv[1]).view(biom.Table)
seqs = pd.Series(t.ids(axis='observation'), index=t.ids(axis='observation'))
seqs_ar = qiime2.Artifact.import_data('FeatureData[Sequence]', seqs)
seqs_ar.save(sys.argv[2])
