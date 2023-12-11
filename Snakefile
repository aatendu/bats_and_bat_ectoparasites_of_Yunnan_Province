"""
File that runs analysis for our bat-bat ectoparasite manuscript.

Written by Alexander Tendu.

Updated 30-Nov-2023.

"""

CONDITIONS = glob_wildcards("reads/{pp_pool}_1.fq.gz").pp_pool
print("pp_pool are: ", CONDITIONS)

dedup_CONDITIONS = glob_wildcards("dedup_reads/{dedup_pool}_1.fq.gz").dedup_pool
print("dedup_pool are: ", dedup_CONDITIONS)

PAIRS = ["g1", "g2", "g3", "g4", "g5", "g6", "g7", "g8", "g9"]
print("g_pairs are: ", PAIRS)

rule all:
   input:
       expand("g_pairs/{g_pair}.co.csv",g_pair=PAIRS)

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
       outdir = "/home/tendu/wanda_cont/trimmed/"
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

#this rule generates a prediction and p-value for the probability of virus seq in each contig. The output is a list of contigs and their scores & p-values
rule dvf_prediction:
    output:
        predictions="dvf_pred/{this_pool}/contigs.fasta_gt1000bp_dvfpred.txt"
    input:
        models="/home/tendu/DeepVirFinder/models/",
        queries=rules.spades_assembly.output.contigs
    params:
        min_query_length=1000,
        outputdir="dvf_pred/{this_pool}/",
        source="/home/tendu/DeepVirFinder/dvf_pred/{this_pool}/contigs.fasta_gt1000bp_dvfpred.txt",
        dest="/home/tendu/wanda_cont/dvf_pred/{this_pool}/contigs.fasta_gt1000bp_dvfpred.txt",
        rece="/home/tendu/wanda_cont/dvf_pred/{this_pool}/"
    threads: 50
    shell:
        """
        cd /home/tendu/DeepVirFinder/
        python dvf.py\
        -i {input.queries}\
        -m {input.models}\
        -o {params.outputdir}\
        -l {params.min_query_length}\
        -c 4
        cp {params.source}  {params.dest}
        cd /home/tendu/wanda_cont/
        """

#this rule removes the dvf predictions with a pvalue greater that .05
rule filter_dvf_output_by_pvalue:
    output:
        record="dvf_filtered_record/{this_pool}.dvf_pred.txt",
        to_use="dvf_filtered_to_use/{this_pool}.dvf_pred.txt"
    input: "dvf_pred/{this_pool}/contigs.fasta_gt1000bp_dvfpred.txt"
    threads: 50
    shell:
        "awk '{{if ($4 < 0.05) {{print$1,$4}} }}' {input} | tee {output.record} | awk '{{print$1}}' > {output.to_use}"

#this rule uses the names of the contigs with p<.05 to subset the original fasta file, retaining the contigs of interest
#the perl function appends a prefix of pool ID to each contig name
rule fetch_fasta_for_dvf_positive:
    output: "dvf_filtered_to_use/{this_pool}_pre_queries.fasta"
    input:
        omnibus_fasta="spades_out/{this_pool}/contigs.fasta",
        names="dvf_filtered_to_use/{this_pool}.dvf_pred.txt"
    shell:
        "seqtk subseq {input.omnibus_fasta} {input.names} | perl -pi -e 's/^>/>{wildcards.this_pool}_/g' > {output}"


#contigs for the pool pairs are coalesced to make a single fasta file.
rule combine_bat_bat_ectoparasite_pool_pairs:
    output:
        g1="g_pairs/g1.fasta",
        g2= "g_pairs/g2.fasta",
        g3= "g_pairs/g3.fasta",
        g4= "g_pairs/g4.fasta",
        g5= "g_pairs/g5.fasta",
        g6= "g_pairs/g6.fasta",
        g7= "g_pairs/g7.fasta",
        g8= "g_pairs/g8.fasta",
        g9= "g_pairs/g9.fasta"

    input:
        b32="dvf_filtered_to_use/b32_pre_queries.fasta",
        b49= "dvf_filtered_to_use/b49_pre_queries.fasta",
        b26= "dvf_filtered_to_use/b26_pre_queries.fasta",
        b28= "dvf_filtered_to_use/b28_pre_queries.fasta",
        b75= "dvf_filtered_to_use/b75_pre_queries.fasta",
        b77= "dvf_filtered_to_use/b77_pre_queries.fasta",
        b30= "dvf_filtered_to_use/b30_pre_queries.fasta",
        b41= "dvf_filtered_to_use/b41_pre_queries.fasta",
        b44= "dvf_filtered_to_use/b44_pre_queries.fasta",
        b116= "dvf_filtered_to_use/b116_pre_queries.fasta",
        b18= "dvf_filtered_to_use/b18_pre_queries.fasta",
        b21= "dvf_filtered_to_use/b21_pre_queries.fasta",
        b66= "dvf_filtered_to_use/b66_pre_queries.fasta",
        b67= "dvf_filtered_to_use/b67_pre_queries.fasta",
        t17= "dvf_filtered_to_use/t17_pre_queries.fasta",
        t18= "dvf_filtered_to_use/t18_pre_queries.fasta",
        t6= "dvf_filtered_to_use/t6_pre_queries.fasta",
        t10= "dvf_filtered_to_use/t10_pre_queries.fasta",
        t1= "dvf_filtered_to_use/t1_pre_queries.fasta",
        t2= "dvf_filtered_to_use/t2_pre_queries.fasta",
        t5= "dvf_filtered_to_use/t5_pre_queries.fasta",
        t13= "dvf_filtered_to_use/t13_pre_queries.fasta",
        t14= "dvf_filtered_to_use/t14_pre_queries.fasta",
        t15= "dvf_filtered_to_use/t15_pre_queries.fasta"

    threads: 16
    shell:
        """
        cat {input.b32} {input.t17} > {output.g1}
        cat {input.b49} {input.t18} > {output.g2}
        cat {input.b26} {input.b28} {input.t6} > {output.g3}
        cat {input.b75} {input.b77} {input.t10} > {output.g4}
        cat {input.b30} {input.t1} > {output.g5}
        cat {input.b41} {input.b44} {input.t2} > {output.g6}
        cat {input.b116} {input.t5} > {output.g7}
        cat {input.b18} {input.b21} {input.t13} > {output.g8}
        cat {input.b66} {input.b67} {input.t14} {input.t15} > {output.g9}

        """

#this rule searches by the blast alg on the ref_viruses_rep_genomes for matches within our combined contig sets g1 to g9
rule blastn_ref_viruses:
    output: "blastn_out/{g_pair}_blastn.csv"
    input: "g_pairs/{g_pair}.fasta"
    params:
        db="/home/tendu/wanda_cont/ref_viruses_may_2023/ref_viruses_rep_genomes",
        evalue=0.00001
    threads: 50
    shell:
        "blastn -db {params.db} -query {input} -out {output} -evalue {params.evalue} -outfmt '10 delim= qseqid qlen sseqid sacc length pident staxid'"

#this rule retrieves the pool prefices to be used later while merged to the blast output. this means the letter, either t or b and puts it in a single column file
rule retrieve_prefix:
    output: "blastn_out/prefix_{g_pair}.txt"
    input: "blastn_out/{g_pair}_blastn.csv"
    shell:
        "awk '{{print substr($0,0,1);}}' {input} > {output}"

#this rule will append the prefix to the blastn output, as the new first column, the rule further converts all tab separations to csv
rule append_prefix:
    output:
        before_sort="blastn_out/before_processing.{g_pair}.csv",
        after_sort="blastn_out/prefixed_{g_pair}_blastn.csv"
    input:
        a="blastn_out/{g_pair}_blastn.csv",
        b="blastn_out/prefix_{g_pair}.txt"
    shell:
        "paste {input.b} {input.a} | sed 's/\t/,/g' | tee {output.before_sort} | sort -k 8,8 -n | awk -F ',' '!seen[$2]++' > {output.after_sort}"

#to see which taxids have matches in both the ectoparasite and bat pools. Use custom python script
rule filter_to_retain_cooccurring_taxa:
    output:
        "g_pairs/{g_pair}.co.csv"
    input:
        "blastn_out/prefixed_{g_pair}_blastn.csv"
    script:
        "scripts/filter_co.py"