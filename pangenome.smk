import pandas as pd
import itertools
import os
import re

configfile: "config/pangenome.yaml"

samples = pd.read_csv(config["samples"], sep='\t')
samples = samples[samples["sample_type"] == "Baseline"]["sample_id"]


genome_derep_ids=pd.read_csv(config["genome_derep_ids"],sep='\t')

checkM=pd.read_csv(config["checkM"],sep='\t')

genome_information = pd.read_csv(config["genome_information"],sep=',')

def get_sample_name(s):
    match = re.search(r'S-\d+', s)
    if match:
        return match.group()
    return None

genome_information = genome_information[ genome_information["genome"].apply(get_sample_name).isin(samples)]

print("number of genomes from relevant samples: "+str(len(genome_information["genome"])))

if config["filter_criteria"] == "combined":
    genome_information["qual"] = genome_information["completeness"]-(5*genome_information["contamination"])
    genome_information_filtered = genome_information[genome_information["qual"] > config["genome_comb_threshold"]]
else:
    genome_information_filtered=genome_information[ (genome_information["completeness"] > config["genome_completeness_threshold"]) &
                            (genome_information["contamination"] < config["genome_contamination_threshold"])]
# if config["filter_criteria"] == "combined":
#     checkM["qual"] = checkM["Completeness"]-(5*checkM["Contamination"])
#     checkM_filtered = checkM[checkM["qual"] > config["genome_comb_threshold"]]
# else:
#     checkM_filtered=checkM[ (checkM["Completeness"] > config["genome_completeness_threshold"]) &
#                             (checkM["Contamination"] < config["genome_contamination_threshold"])]

genome_list = genome_derep_ids[ genome_derep_ids["genome"].isin(genome_information_filtered["genome"].str.replace(".fasta",""))]

print("number of genomes after filtering by qual: "+str(len(genome_list["genome"])))

gtdb_bac = pd.read_csv(config["gtdb_bac"], sep='\t')
gtdb_ar = pd.read_csv(config["gtdb_ar"], sep='\t')

gtdb_dat = pd.concat([gtdb_bac[["user_genome","classification"]],
                        gtdb_ar[["user_genome","classification"]]], ignore_index=True)

def extract_species_name(s):
    match = re.search(r's__(.*)', s)
    if match:
        return match.group(1).replace(" ", "_")
    return ''

gtdb_dat["species"] = gtdb_dat["classification"].apply(extract_species_name)
# gtdb_dat.drop("classification", axis=1, inplace=True)

genome_list = pd.merge(genome_list, gtdb_dat, left_on="MAG", right_on="user_genome", how='left')
genome_list.drop("user_genome", axis=1, inplace=True)

genome_list.loc[genome_list["species"] == "", "species"] = genome_list["MAG"]

print(genome_list['species'].value_counts().sort_values(ascending=False))

## Analyze Phascolarctobacterium A succinatutens A together with the other ones.
genome_list["species"] = genome_list["species"].replace("Phascolarctobacterium_A_succinatutens_A", "Phascolarctobacterium_A_succinatutens")

## Restrict the analysis to these species
genomes=[ "Phascolarctobacterium_A_sp900544885", "Phascolarctobacterium_A_succinatutens", "Phascolarctobacterium_A_sp900552855",
           "Phascolarctobacterium_A_sp900541915", "Phascolarctobacterium_A_sp900552005", "Phascolarctobacterium_faecium",
           "Phascolarctobacterium_sp900545535", "Phascolarctobacterium_sp000436095"]

genome_list = genome_list[genome_list["species"].isin(genomes)]
# raise Exception("Enough")

combined_species = genome_list.copy()
combined_species["species"] = "Phascolarcto"

print("number of genomes before adding combined: "+str(len(genome_list["genome"])))

genome_list = pd.concat([genome_list, combined_species], ignore_index=True)

genome_list = genome_list.groupby('MAG').filter(lambda x: len(x) >= config["min_genomes"])

print("number of genomes after filtering by freq: "+str(len(genome_list["genome"])))


rule all:
    input:
        # [ os.path.join(config["wd"], "pangenome", "prokka", row["species"], row["genome"]+".gff") for _, row in genome_list.iterrows()],
        expand(os.path.join(config["wd"],"pangenome","roary", "{species}", "summary_statistics.txt"), species=set(genome_list["species"])),
        expand(os.path.join(config["wd"],"pangenome","plots","{species}","pangenome_frequency.png"), species = set(genome_list["species"])),
        expand(os.path.join(config["wd"],"pangenome", "roary", "{species}", "gene_presence_long.tsv"), species = set(genome_list["species"]))

rule annotate_genomes:
    input:
        [ os.path.join(config["wd"], "pangenome", "dram", row["species"], row["genome"], "annotations.tsv") for _, row in genome_list.iterrows()]

def get_genome_path(wildcard):
    sample=re.search(r'(.+?)(?=_metabat|_maxbin)', wildcard.genome)
    sample=sample.group(1)
    return os.path.join(config["atlas"],sample,"binning","DASTool", "bins", wildcard.genome+".fasta")

# def get_genus(species):
#     classification = gtdb_dat[gtdb_dat["species"] == species]["classification"][0].split(";")
#     genus = [ this_one for this_one in classification if this_one[:3] == "g__"][0].str.replace("g__","")
#     return genus

rule prokka:
    input:
        genome = lambda w: get_genome_path(w)
    output:
        prokka_res = os.path.join(config["wd"], "pangenome", "prokka", "{species}", "{genome}.gff"),
    params:
        outdir = os.path.join(config["wd"],"pangenome","prokka","{species}",""),
        prefix = lambda w: f'{w.genome}',
        min_contig_length = 200,
        # genus = lambda w: get_genus(w.species),
        log_file = os.path.join(config["wd"], "pangenome", "prokka", "{species}", "{genome}.log")
    log:
        os.path.join(config["wd"],"logs","prokka","{species}_{genome}.log")
    conda:
        "envs/prokka.yaml"
    shadow: "minimal"
    threads:
        8
    shell:
        """
            prokka {input.genome} \
                --outdir {params.outdir} \
                --prefix {params.prefix} \
                --gcode 11 \
                --cpus {threads} \
                --mincontiglen {params.min_contig_length} \
                --force > /dev/null 2>&1 
            mv {params.log_file} {log}
        """
                # --genus {wildcards.genus} \
                # --usegenus \

rule flag_prokka:
    input:
        prokka_output=[ os.path.join(config["wd"], "pangenome", "prokka", row["species"], row["genome"]+".gff") for _, row in genome_list.iterrows()]
    output:
        flag=temp(os.path.join(config["wd"], "pangenome", "prokka_flag"))
    shell:
        "touch {output.flag}"

rule dram:
    input:
        genome = lambda w: get_genome_path(w)
    output:
        dram_annotations = os.path.join(config["wd"], "pangenome", "dram", "{species}", "{genome}","annotations.tsv"),
        dram_gff = os.path.join(config["wd"], "pangenome", "dram", "{species}", "{genome}", "genes.gff"),
        dram_fasta = os.path.join(config["wd"], "pangenome", "dram", "{species}", "{genome}", "scaffolds.fna")
    params:
        out_dir=os.path.join(config["wd"], "pangenome", "dram", "{species}", "{genome}", "")
    log:
        os.path.join(config["wd"],"logs","dram","{species}_{genome}.log"),
        # os.path.join()"logs/{stringency}/dramv/dramv.log"
    conda:
        "envs/dram.yaml"
    threads:
        1
    benchmark:
        os.path.join(config["wd"], "benchmarks", "dram", "{species}", "{genome}.txt")
    shell:
        """
        
        [ -d {params.out_dir} ]  && rm {params.out_dir} -r

        DRAM.py annotate \
            --input_fasta {input.genome} \
            --output_dir {params.out_dir} \
            --use_uniref \
            --threads {threads} \
            --verbose \
            >& {log}
        """

rule get_dram_gffs:
    input:
        in_gff=rules.dram.output.dram_gff,
        in_fasta=rules.dram.output.dram_fasta,
    output:
        out_gff=os.path.join(config["wd"], "pangenome","dram","{species}","gffs","{genome}.gff"),
    conda:
        "envs/seqtk.yaml"
    shell:
        """
            echo "##gff-version 3" > {output.out_gff}
            grep -v "##gff-version 3" {input.in_gff} >> {output.out_gff}
            echo "##FASTA" >> {output.out_gff}
            seqtk seq -C {input.in_fasta} >> {output.out_gff}
        """

rule flag_dram:
    input:
        dram_output=[ os.path.join(config["wd"], "pangenome", "dram", row["species"], row["genome"], "genes.gff") for _, row in genome_list.iterrows()]
    output:
        flag=temp(os.path.join(config["wd"], "pangenome", "dram_flag"))
    shell:
        "touch {output.flag}"

def get_species_gffs(wildcard, annotation_source):
    if annotation_source == "prokka":
        gff_paths = [ os.path.join(config["wd"],"pangenome","prokka",wildcard.species,genome+".gff") 
                    for genome in genome_list[genome_list["species"] == wildcard.species]["genome"]]
    elif annotation_source == "dram":
        gff_paths = [ os.path.join(config["wd"],"pangenome","dram",wildcard.species,"gffs",genome+".gff") 
                    for genome in genome_list[genome_list["species"] == wildcard.species]["genome"]]
    return gff_paths

annotation_type = "dram"

rule roary:
    input:
        annotation_flag = os.path.join(config["wd"], "pangenome", annotation_type+"_flag"),
        gffs = lambda w: get_species_gffs(w, annotation_type)
    output:
        summary_stats=os.path.join(config["wd"],"pangenome", "roary", "{species}", "summary_statistics.txt"),
        gene_presence=os.path.join(config["wd"],"pangenome", "roary", "{species}", "gene_presence_absence.csv"),
        alignment=os.path.join(config["wd"],"pangenome", "roary", "{species}", "core_gene_alignment.aln"),
        clustered_proteins=os.path.join(config["wd"],"pangenome", "roary", "{species}", "clustered_proteins"),
        # alignment=config["wd"]+"pangenome/{tax_db}/{genus}/roary/core_gene_alignment.aln"
    params:
        tmp_folder=os.path.join(config["wd"],"pangenome", "roary", "{species}","roary_tmp",""),
        out_folder=os.path.join(config["wd"],"pangenome", "roary", "{species}",""),
    log:
        os.path.join(config["wd"],"logs","roary","{species}.log")
    conda:
        "envs/roary.yaml"
    threads:
        8
    shell:
        """
            roary -p {threads} \
                -f {params.tmp_folder} \
                -e \
                --mafft \
                -i 75 \
                -t 11 \
                {input.gffs} &> {log}
            
            mv {params.tmp_folder}* {params.out_folder}.
            rm {params.tmp_folder} -r
        """

rule core_genome_tree:
    input:
        alignment=rules.roary.output.alignment
    output:
        core_genome_tree=os.path.join(config["wd"],"pangenome","roary","{species}","core_gene_alignment.aln.iqtree"),
        tree=os.path.join(config["wd"],"pangenome","roary","{species}","core_gene_alignment.aln.treefile"),
    log:
        os.path.join(config["wd"],"logs","roary","{species}.log")
    conda:
        "envs/iqtree.yaml"
    shell:
        """
            iqtree -s {input.alignment} -m GTR+F+R7 -redo &> {log}
        """

rule annotate_roary_long:
    input:
        wide_roary=rules.roary.output.clustered_proteins,
        annotations=lambda w: [ os.path.join(config["wd"], "pangenome", "dram", w.species, genome,"annotations.tsv") 
                                for genome in genome_list[genome_list["species"] == w.species]["genome"] ]
    output:
        annotated_groups=os.path.join(config["wd"],"pangenome", "roary", "{species}", "gene_presence_long.tsv"),
    run:
        annotations = pd.concat([ pd.read_csv(file, delimiter="\t", index_col=0) for file in input.annotations ])
        annotations.reset_index(inplace=True)
        annotations.rename(columns={'index': 'protein_id'}, inplace=True)
        
        sep_func = lambda x: pd.Series([i for i in x.split(': ') + [''] if i])

        groups = pd.read_csv(input.wide_roary, delimiter=': |\t', header=None, engine='python')
        groups.rename(columns={0: 'Identifier'}, inplace=True)
        groups['protein_id'] = groups.iloc[:, 1:].apply(lambda x: list(x.dropna()), axis=1)
        # groups.rename(columns={0: 'DB', 1: 'Identifier'}, inplace=True)
        # groups['protein_id'] = groups.iloc[:, 2:].apply(lambda x: list(x.dropna()), axis=1)
        groups_long = groups[["Identifier", "protein_id"]].explode('protein_id')
        print(groups_long)

        annotations_tab = pd.merge(groups_long, 
                                    annotations, 
                                    left_on='protein_id', 
                                    right_on='protein_id', 
                                    how='left')
        
        annotations_tab.to_csv(output.annotated_groups, sep="\t")
