"""
First file that runs our analysis for the bat-bat ectoparasite study, in the revised version of the manuscript
(from Assembled contigs onwards)

written by A.A.Tendu

Updated 20 May 2024
"""

CONDITION = glob_wildcards("filtered_reads/wanted_pp_reads_{this_pool}_fwd.fq.gz").this_pool
#print("this_pool are: ", CONDITION) #if not commented out, snakemake cannot execute the --rulegraph command that produces a flowchart

B_T_PAIRS = ["g1b","g1t","g2b", "g2t", "g3b", "g3t", "g4b", "g4t", "g5b", "g5t", "g6b", "g6t", "g7b", "g7t", "g8b", "g8t", "g9b", "g9t"]

rule all:
   input:
       first= expand("merged_cpm_blastn/merged_{this_pool}.csv",this_pool=CONDITION),
       second = expand("new_pools/{this_g}.csv",this_g=B_T_PAIRS)

#this rule filters the master fasta file to retain the sequences that are longer than 500bp
rule filter_by_size:
    output: "over_500/{this_pool}_a500.fasta"
    input: "spades_out/{this_pool}.fasta"
    shell:
        """awk -v n=500 '/^>/{{ if(l>n) print b; b=$0;l=0;next }} {{l+=length;b=b ORS $0}}END{{if(l>n) print b }}' {input} > {output}"""

#cheromydarling
#this rule appends a prefix to the start of each sequence according to the pool IDs
rule rename_the_sequences:
    output: "renamed/{this_pool}.fasta"
    input: "over_500/{this_pool}_a500.fasta"
    shell:
        "cat {input} | perl -pi -e 's/^>/>{wildcards.this_pool}_/g' > {output}"

#this rule builds a bowtie2 index using all contigs longer than 500bp for each pool
rule build_bowtie_index:
    output:
        one = "btindex/{this_pool}.1.bt2",
        two = "btindex/{this_pool}.2.bt2",
        three = "btindex/{this_pool}.3.bt2",
        four = "btindex/{this_pool}.4.bt2",
        five = "btindex/{this_pool}.rev.1.bt2",
        six = "btindex/{this_pool}.rev.2.bt2"
    input:
        fasta =rules.rename_the_sequences.output
    threads: 24
    shell:
        "bowtie2-build -f {input.fasta} btindex/{wildcards.this_pool}"

rule map_the_reads:
    output:
        sam = "{this_pool}.sam"
    input:
        reads_f = "filtered_reads/wanted_pp_reads_{this_pool}_fwd.fq.gz",
        reads_r = "filtered_reads/wanted_pp_reads_{this_pool}_rev.fq.gz",
        index =rules.build_bowtie_index.output
    params:
        bow2_index = "btindex/{this_pool}"
    threads: 24
    shell:
        """
        bowtie2\
        --sensitive-local\
        -x {params.bow2_index}\
        -2 {input.reads_r}\
        -1 {input.reads_f}\
        -t\
        -S {output.sam}
        """
rule convert_sam_to_bam:
    output: "{this_pool}.bam"
    input:rules.map_the_reads.output.sam
    threads: 24
    shell:
        "samtools view -b -o {output} {input}"

rule sort_bam_file:
    output: "{this_pool}_sorted.bam"
    input:rules.convert_sam_to_bam.output
    threads: 24
    shell:
        "samtools sort -o {output} {input}"

rule obtain_read_stats:
    output: "coverage/{this_pool}_read_stats.tsv"
    input:rules.sort_bam_file.output
    threads: 24
    shell:
        "samtools coverage -o {output} {input}"

rule run_blastn_on_all_contigs:
    output: "blastn_out/raw_{this_pool}.csv"
    input: "renamed/{this_pool}.fasta"
    params:
        db = "/home/tendu/new_pplus/virus_db/ref_viruses_rep_genomes"
    shell:
        "blastn -db {params.db} -query {input} -out {output} -evalue 0.001 -outfmt '10 delim= qseqid qlen staxid pident length sacc' -subject_besthit -max_target_seqs 5"

# Here duplicate hits refers to a single contig showing more than one hit. we retain the first
rule remove_duplicate_hits:
    output:
         clean = "blastn_out/clean/clean_{this_pool}.csv",
         lineage = "lineage/{this_pool}_lineage.ssv"
    input: "blastn_out/raw_{this_pool}.csv"
    shell:
        "uniq -w 13 {input} | tee {output.clean} | awk -F ',' '{{print $3}}' | taxonkit lineage --delimiter ',' > {output.lineage}"

##next step is to unpair the g_pairs and run the blastn on all pools separately. then remove duplicates (kmbienes_script). then add a cpm generation
## rule. add the cpm column to the blastn_outputs. count the richness, abundance, of predicted viruses ets, for g_pairs

#this rule creates a column: cpm_reads from the numreads column produced by samtools coverage (cpm=counts per million) i.e. reads*10^6/total_reads_in_pool
rule calculate_cpm_for_mapped_reads:
    output: "normalised_coverage/cpm_{this_pool}.csv"
    input: "coverage/{this_pool}_read_stats.tsv"
    script:
        "scripts/normalize.py"
#this rule will merge the two dataframes, one containing the blastn output and the other containing cpm (normalised reads)
rule merge_blastn_with_normalised_reads:
    output: "merged_cpm_blastn/merged_{this_pool}.csv"
    input:
        cpm_coverage = "normalised_coverage/cpm_{this_pool}.csv",
        blastn_clean = "blastn_out/clean/clean_{this_pool}.csv"
    script:
        "scripts/add_cpm_to_blastn.py"

#this rule will retain the hit the query with the highest cpm_reads for each staxid in all pools separately
rule remove_duplicate_queries_within_taxid:
    output: "no_duplicate_staxid/merged_dedup_{this_pool}.csv"
    input: "merged_cpm_blastn/merged_{this_pool}.csv"
    threads: 16
    script:
        "scripts/remove_dupl_taxids.py"

#this rule should combine those pools that will be regarded as a single pool in subsequent analysis
rule combine_new_g_pools:
    output:
        g3b = "pre_pre_new_pools/g3b.csv",
        g4b = "pre_pre_new_pools/g4b.csv",
        g6b = "pre_pre_new_pools/g6b.csv",
        g8b = "pre_pre_new_pools/g8b.csv",
        g9b = "pre_pre_new_pools/g9b.csv",
        g9t = "pre_pre_new_pools/g9t.csv"
    input:
        b26 = "no_duplicate_staxid/merged_dedup_b26.csv",
        b28 = "no_duplicate_staxid/merged_dedup_b28.csv",
        b75 = "no_duplicate_staxid/merged_dedup_b75.csv",
        b77 = "no_duplicate_staxid/merged_dedup_b77.csv",
        b41 = "no_duplicate_staxid/merged_dedup_b41.csv",
        b44 = "no_duplicate_staxid/merged_dedup_b44.csv",
        b18 = "no_duplicate_staxid/merged_dedup_b18.csv",
        b21 = "no_duplicate_staxid/merged_dedup_b21.csv",
        b66 = "no_duplicate_staxid/merged_dedup_b66.csv",
        b67 = "no_duplicate_staxid/merged_dedup_b67.csv",
        t14 = "no_duplicate_staxid/merged_dedup_t14.csv",
        t15 = "no_duplicate_staxid/merged_dedup_t15.csv"
    threads: 16
    shell:
        """
        cat {input.b26} {input.b28} > {output.g3b}
        cat {input.b75} {input.b77} > {output.g4b}
        cat {input.b41} {input.b44} > {output.g6b}
        cat {input.b18} {input.b21} > {output.g8b}
        cat {input.b66} {input.b67} > {output.g9b}
        cat {input.t14} {input.t15} > {output.g9t}
        """
 #to remove the additional header line resulting from the cat command
rule remove_extra_header_line:
    output : "pre_new_pools/{this_g}.csv"
    input : "pre_pre_new_pools/{this_g}.csv"
    threads: 16
    shell:
        "awk '!a[$0]++' {input} > {output}"
#some pools dont need merging
rule bring_the_singlet_pools_by_renaming:
    output:
        g1b = "new_pools/g1b.csv",
        g1t = "new_pools/g1t.csv",
        g2b = "new_pools/g2b.csv",
        g2t = "new_pools/g2t.csv",
        g3t = "new_pools/g3t.csv",
        g4t = "new_pools/g4t.csv",
        g5b = "new_pools/g5b.csv",
        g5t = "new_pools/g5t.csv",
        g6t = "new_pools/g6t.csv",
        g7b = "new_pools/g7b.csv",
        g7t = "new_pools/g7t.csv",
        g8t = "new_pools/g8t.csv"
    input:
        b32 = "no_duplicate_staxid/merged_dedup_b32.csv",
        t17 = "no_duplicate_staxid/merged_dedup_t17.csv",
        b49 = "no_duplicate_staxid/merged_dedup_b49.csv",
        t18 = "no_duplicate_staxid/merged_dedup_t18.csv",
        t6 = "no_duplicate_staxid/merged_dedup_t6.csv",
        t10 = "no_duplicate_staxid/merged_dedup_t10.csv",
        b30 = "no_duplicate_staxid/merged_dedup_b30.csv",
        t1 = "no_duplicate_staxid/merged_dedup_t1.csv",
        t2 = "no_duplicate_staxid/merged_dedup_t2.csv",
        b116 = "no_duplicate_staxid/merged_dedup_b116.csv",
        t5 = "no_duplicate_staxid/merged_dedup_t5.csv",
        t13 = "no_duplicate_staxid/merged_dedup_t13.csv"
    shell:
        """
        mv {input.b32} {output.g1b}
        mv {input.t17} {output.g1t}
        mv {input.b49} {output.g2b}
        mv {input.t18} {output.g2t}
        mv {input.t6} {output.g3t}
        mv {input.t10} {output.g4t}
        mv {input.b30} {output.g5b}
        mv {input.t1} {output.g5t}
        mv {input.t2} {output.g6t}
        mv {input.b116} {output.g7b}
        mv {input.t5} {output.g7t}
        mv {input.t13} {output.g8t}
        """
#also remove duplicate taxids from the new g_pools
rule remove_dupl_in_new_g_groups:
    output: "new_pools/{this_sub_g}.csv"
    input: "pre_new_pools/{this_sub_g}.csv"
    threads: 12
    script:
        "scripts/remove_dupl_taxids.py"

