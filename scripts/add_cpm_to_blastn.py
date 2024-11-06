## this script brings the cpm_reads column into the blastn output

import pandas as pd

cpm_coverage = pd.read_csv(snakemake.input[0])

cpm_coverage = cpm_coverage.rename(columns={'#rname': 'qseqid'})

blastn_clean = pd.read_csv(snakemake.input[1], names=["qseqid", "qlen", "staxid", "pident", "length", "sacc"])

new_coverage = cpm_coverage.drop(['startpos', 'endpos', 'covbases', 'coverage', 'meandepth', 'meanbaseq', 'meanmapq'], axis=1)

new_coverage = new_coverage.merge(blastn_clean, left_on='qseqid', right_on='qseqid')

new_coverage.to_csv(snakemake.output[0], index=False)