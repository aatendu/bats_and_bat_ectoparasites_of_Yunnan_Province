"""meant to see which taxids have matches in both the bat and the ectoparasite pools"""
import pandas as pd

mydf = pd.read_csv(snakemake.input[0], names=["pool", "qseqid", "qlen", "sseqid", "sacc", "length", "pident", "staxid"])

grps = mydf.groupby(['staxid'])['pool'].apply(set)

filtered = grps[grps.map(len) > 1].index

res = mydf[mydf.set_index(['staxid']).index.isin(filtered)]

res.to_csv(snakemake.output[0], index=False)
