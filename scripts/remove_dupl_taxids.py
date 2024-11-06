#this script should retain only one row for all row sets sharing one taxid
import pandas as pd

master = pd.read_csv(snakemake.input[0])

master_new = master.sort_values('cpm_reads', ascending=False).drop_duplicates('staxid').sort_index()

master_new.to_csv(snakemake.output[0], index=False)