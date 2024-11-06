##this script only adds the counts per million column

import pandas as pd

df = pd.read_csv(snakemake.input[0], sep='\t')

df['cpm_reads'] = (df['numreads']*1000000)/ df['numreads'].sum()

df.to_csv(snakemake.output[0], index=False)