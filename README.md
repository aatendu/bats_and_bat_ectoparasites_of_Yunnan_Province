## Viromes of Obligate Parasites and their Hosts: The Case of Bats and Bat Ectoparasites

### Authors
Alexander Tendu<sup>1,2,3</sup>, Ruiya Li<sup>1,3</sup>, Yakhouba Kane<sup>1,3</sup>, Betty Nalikka<sup>1,2,3</sup>, Victor Omondi<sup>2</sup>, Kathrina Mae Bienes<sup>2</sup>, Nicolas Berthet<sup>4,5</sup>, and Gary Wong<sup>1,2,#</sup>.

<sup>1</sup>Institut Pasteur of Shanghai-CAS; <sup>2</sup>Institut Pasteur du Cambodge, Cambodia; <sup>3</sup>University of Chinese Academy of Sciences; <sup>4</sup>Institut Pasteur - EPVO - Epidémiologie et Physiopathologie des Virus Oncogenes, France; <sup>5</sup>Institut Pasteur, Unité Environnement et Risque Infectieux, Cellule d’Intervention Biologique d’Urgence, France.

### Abstract

Blood feeding ectoparasites of bats have been found to contain insect-specific and vertebrate-infecting viruses, some of which are relatives of human infecting viruses. While some of these viruses have been speculated as being of bat origin, only few have been shown to co-occur in their specific bat hosts. To investigate the presence/absence of co-occurring viruses between these species, we compared virus sequences from bats and their blood feeding ectoparasites to investigate the presence/absence of co-occurring viruses between these species. Although vertebrate-infecting viruses were not found to co-occur between both groups, our analysis showed that shared bacteriophages may illuminate genomic and non-genomic characteristics influencing the chance of co-occurrence in these two systems. We observed that the abundance for a majority of co-occurring viruses was higher in bat ectoparasites than in bats. We found that genome length was an important predictor of abundance in either bats or their ectoparasites, whereby an increase in genome length was associated with a greater abundance of co-occurring viruses in bats rather than in their ectoparasites. Our findings illuminate factors that may influence the outcome of vector competence investigation in bat ectoparasites

### About this repository
This repository contains a workflow that should reproduce the data in the manuscript. The complete workflow is contained in [Snakefile](Snakefile).

In brief, the workflow consists of these steps;
1. Fwd & Rev reads are trimmed for adaptor removal and quality threshold filtration.
2. Host-associated sequences are removed by filtering against the default SILVA db (smr_v4.3_default_db.fasta).
3. Reads are assembled into contigs using SPAdes.
4. DeepVirFinder produces a prediction (score and p-value) indicating the probability that each contig is of viral origin.
5. Contigs from ectoparasites and their hosts (bats) with p-values less than .05 are coalesced into a single file that is later used as query for this particular bat-bat ectoparasite pool pair.
6. The RefSeq Viruses rep genomes db is used with blastn to produce matches to each query file (nine in total)
7. Subject sequences in the RefSeq Viruses rep genomes db that have matches in both the ectoparasite and its corresponding bat pool are retained for further analysis.

The workflow may be run in this [conda environment](environment.yaml).
