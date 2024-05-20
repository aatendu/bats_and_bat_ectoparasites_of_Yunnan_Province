"""
First file that runs our analysis for the bat-bat ectoparasite study, in the revised version of the manuscript
(from read processing to Assembly)

written by A.A.Tendu

Updated 20 May 2024
"""

CONDITIONS = glob_wildcards("reads/{pp_pool}_1.fq.gz").pp_pool
print("pp_pool are: ", CONDITIONS)

dedup_CONDITIONS = glob_wildcards("dedup_reads/{dedup_pool}_1.fq.gz").dedup_pool
print("dedup_pool are: ", dedup_CONDITIONS)

this_pool_CONDITIONS = CONDITIONS
print("this_pool are: ", this_pool_CONDITIONS)

rule all:
   input:
       expand("spades_out/{this_pool}/contigs.fasta",this_pool=this_pool_CONDITIONS)

#read trimming by trim galore
rule trim_reads:
   output:
       r1 = "trimmed/{pp_pool}_1_val_1.fq.gz",
       r2 = "trimmed/{pp_pool}_2_val_2.fq.gz"
   input:
       r1 = "reads/{pp_pool}_1.fq.gz",
       r2 = "reads/{pp_pool}_2.fq.gz"
   params:
       q = 29,
       stringency = 3,
       length = 100,
       max_n = 50,
       outdir = "/home/tendu/new_pplus/trimmed/"
   threads: 50
   shell:
       """
         trim_galore --illumina\
           -o {params.outdir}\
           -q {params.q}\
           --length {params.length}\
           --stringency {params.stringency}\
           --max_n {params.max_n}\
           --paired {input.r1} {input.r2}
       """

#read trimming by fastp
rule fastp_trimming:
    output:
        r1 = "trimmed/{dedup_pool}_1_val_1.fq.gz",
        r2 = "trimmed/{dedup_pool}_2_val_2.fq.gz"
    input:
        r1 = "dedup_reads/{dedup_pool}_1.fq.gz",
        r2 = "dedup_reads/{dedup_pool}_2.fq.gz"
    params:
        qual_thresh = 29,
        length_thresh = 100,
        dup_accuracy = 16
    threads: 50
    shell:
        """
        fastp\
        --in1 {input.r1}\
        --in2 {input.r2}\
        --out1 {output.r1}\
        --out2 {output.r2}\
        --detect_adapter_for_pe\
        --dedup_calc_accuracy {params.dup_accuracy}\
        -q {params.qual_thresh}\
        -l {params.length_thresh}
        """

#this rule will filter the reads to remove host related sequences by searching against the SILVA DB
rule silva_filtration:
    output:
        wantedf = "silva_filtered/wanted_{this_pool}_fwd.fq.gz",
        wantedr = "silva_filtered/wanted_{this_pool}_rev.fq.gz",
        unwantedf = "silva_unwanted/unwanted_{this_pool}_fwd.fq.gz",
        unwantedr = "silva_unwanted/unwanted_{this_pool}_rev.fq.gz"
    input:
        ref1="silva_db/smr_v4.3_default_db.fasta",
        reads1 ="trimmed/{this_pool}_1_val_1.fq.gz",
        reads2 = "trimmed/{this_pool}_2_val_2.fq.gz"
    threads: 50
    params:
        ali = "silva_unwanted/unwanted_{this_pool}",
        filt = "silva_filtered/wanted_{this_pool}",
        kvdb = "silva_filtered/{this_pool}/"
    shell:
        """
        sortmerna\
        --ref {input.ref1}\
        --reads {input.reads1}\
        --reads {input.reads2}\
        --aligned {params.ali}\
        --other {params.filt}\
        --kvdb {params.kvdb}\
        --fastx\
        --paired_out\
        --out2
        """

#assembly using spades;
rule spades_assembly:
    output:
        contigs = "spades_out/{this_pool}/contigs.fasta",
        scaffolds = "spades_out/{this_pool}/scaffolds.fasta"
    input:
        right1 ="silva_filtered/wanted_{this_pool}_fwd.fq.gz",
        left1 ="silva_filtered/wanted_{this_pool}_rev.fq.gz"
    params:
        outputdir = "spades_out/{this_pool}/"
    threads: 50
    shell:
        """
        spades.py -1 {input.right1}\
        -2 {input.left1}\
        -o {params.outputdir}\
        -k 21,33,55\
        --careful\
        --disable-gzip-output
        """
