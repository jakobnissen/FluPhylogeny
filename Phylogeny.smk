import os

SNAKEDIR = os.path.dirname(workflow.snakefile)
JULIA_COMMAND = f"JULIA_LOAD_PATH='{SNAKEDIR}' julia --startup-file=no"

def list_dir(path):
    return list(filter(lambda x: x != ".DS_Store", os.listdir(path)))

###############
# Parse config
###############
KNOWN_CONFIGS = {'ref', 'host', 'consensus', 'phylocons'}
for key in config:
    if key not in KNOWN_CONFIGS:
        raise KeyError(
            f"Config \"{key}\" is not a known Consensus pipeline config. "
            "Check spelling. Known configs are: "
            f"{','.join(KNOWN_CONFIGS)}"
        )

# Appalingly, Snakemake automatically converts input strings to int/float
# if they are parseable as such. Convert them back to strings
for (k, v) in config.items():
    config[k] = str(v)

# Get refdir
if "ref" not in config:
    raise KeyError("You must supply reference directory: '--config ref=/path/to/ref'")

TOP_REF_DIR = os.path.abspath(os.path.expanduser(config["ref"]))
if not os.path.isdir(TOP_REF_DIR):
    raise NotADirectoryError(TOP_REF_DIR)

# Get host
if "host" not in config:
    raise KeyError("You must supply host: '--config host=human'")
HOST = config["host"]

# Check refdir has a host subdir
possible_hosts = set(list_dir(TOP_REF_DIR))
if HOST not in possible_hosts:
    raise KeyError(f"Directory for host {HOST} not found in {TOP_REF_DIR}")

# We expect human viruses to be quite close to the reference viruses, and any virus
# significantly different from the seasonal references is of interest, so
# we set a high thresholds.
# Empirically, fairly distant subtypes can be 87% identical on NA.
MINIMUM_IDENTITY = 0.96 if HOST == 'human' else 0.8

REFDIR = os.path.join(TOP_REF_DIR, HOST)
REFOUTDIR = os.path.join(os.path.dirname(TOP_REF_DIR), "refout", "phylo", HOST)

# Get consensus dir
PHYLOCONS = "phylocons" in config
if not (PHYLOCONS ^ ("consensus" in config)):
    raise KeyError(
        "You must supply either \"consensus\" config or \"phylocons\" config: "
        "'--config consensus=/path/to/consensus'"
    )

if PHYLOCONS:
    CONSENSUS_DIR = os.path.abspath(os.path.expanduser(config["phylocons"]))
else:
    CONSENSUS_DIR = os.path.abspath(os.path.expanduser(config["consensus"]))

if not os.path.isdir(CONSENSUS_DIR):
    raise NotADirectoryError(CONSENSUS_DIR)

ALL_SEGMENTS = sorted(list_dir(os.path.join(REFDIR, "segments")))

###########################
# This is the function that triggers the checkpoint. By accessing the blastout
# checkpoint, the DAG "short-circuits" all the middle rules until the file
# has been completed.
def all_inputs(wildcards):
    # This is the function that gets all the segment/type combinations into the DAG
    trigger = checkpoints.genotypes.get()
    files = ["genotypes.txt"]
    groups = [p[:-4].partition('_') for p in list_dir("tmp/catgroups")]
    files.extend([f"trees/{s}/{g}.pdf" for (s, _, g) in groups])
    return files

rule all:
    input: all_inputs
    output: "phylogeny_versions.txt"
    params: SNAKEDIR
    shell:
        "cp {params:q}/copy_readme.md README_PHYLOGENY.md && "
        "(git -C {params:q} rev-parse --short HEAD > {output} || true) && "
        "julia -v >> {output}"

###########################
# Before checkpoint: BLASTn
###########################
rule instantiate:
    output: touch(REFOUTDIR +  "/phylo_instantiated")
    params: JULIA_COMMAND
    shell: "{params} -e 'using Pkg; Pkg.resolve(); Pkg.instantiate()'"

# Cat all reference segments together, and all tree groups
rule cat_ref:
    output:
        marker=REFOUTDIR + "/phylo_ref_cat_done",
        cat=expand(REFOUTDIR + "/segments/{segment}.fna", segment=ALL_SEGMENTS)
    params:
        juliacmd=JULIA_COMMAND,
        scriptpath=f"{SNAKEDIR}/scripts/cat_ref.jl",
        refdir=REFDIR,
        known_genotypes=os.path.join(REFDIR, "genotypes.tsv"),
        refoutdir=REFOUTDIR
    shell: 
        "{params.juliacmd} {params.scriptpath:q} "
        "{params.refdir:q} {params.known_genotypes:q} {params.refoutdir:q}"

rule makeblastdb:
    input: REFOUTDIR + "/segments/{segment}.fna"
    output:
        nin=REFOUTDIR + "/segments/{segment}.fna" + ".nin",
        nhr=REFOUTDIR + "/segments/{segment}.fna" + ".nhr",
        nsq=REFOUTDIR + "/segments/{segment}.fna" + ".nsq"
    log: "tmp/log/makeblastdb/{segment}.log"
    shell: "makeblastdb -in {input:q} -dbtype nucl 2> {log}"

rule gather_cons:
    input: REFOUTDIR +  "/phylo_instantiated"
    output: temp(expand("tmp/catcons/{segment}.fna", segment=ALL_SEGMENTS))
    params:
        juliacmd=JULIA_COMMAND,
        scriptpath=f"{SNAKEDIR}/scripts/gather_consensus.jl",
        segment_dir=os.path.join(REFDIR, "segments"),
        consensus_dir=CONSENSUS_DIR,
        phylocons=str(PHYLOCONS).lower()
    shell: "{params.juliacmd} {params.scriptpath:q} tmp/catcons {params.segment_dir:q} {params.consensus_dir:q} {params.phylocons}"

rule blastn:
    input:
        nin=REFOUTDIR + "/segments/{segment}.fna" + ".nin",
        nhr=REFOUTDIR + "/segments/{segment}.fna" + ".nhr",
        nsq=REFOUTDIR + "/segments/{segment}.fna" + ".nsq",
        fna="tmp/catcons/{segment}.fna"
    output: "tmp/blast/{segment}.blastout"
    params:
        basename=REFOUTDIR + "/segments/{segment}.fna"
    shell: "blastn -query {input.fna} -db {params.basename} -outfmt '6 qacc sacc qcovhsp pident bitscore' > {output}"

# Purpose: Determine the subtype for each segment in ALL_SEGMENTS.
# write the output to tmp/clades.txt
# write human-readable report to genotypes.txt
checkpoint genotypes:
    input:
        blast=expand("tmp/blast/{segment}.blastout", segment=ALL_SEGMENTS),
        cons=expand("tmp/catcons/{segment}.fna", segment=ALL_SEGMENTS)
    output:
        genotypes="genotypes.txt"
        # Also tmp/cat/{segment}_{clade}.fna for every seg/type found
        # but the pairs ar unknown at this point, so it's not part of the rule
    params:
        juliacmd=JULIA_COMMAND,
        scriptpath=f"{SNAKEDIR}/scripts/parse_blast.jl",
        catconsdir="tmp/catcons", # corresponds to input.cons
        outconsdir="tmp/catgroups", # dir of output .fna files
        jlspath="tmp/genotypes.jls.gz", # path of internal genotypes object
        inconsdir=CONSENSUS_DIR, # we use this to load a list of sample names
        blastdir="tmp/blast", # corresponds to input.blast
        tree_groups=os.path.join(REFDIR, "tree_groups.txt"),
        known_genotypes=os.path.join(REFDIR, "genotypes.tsv"),
        minid=MINIMUM_IDENTITY,
        phylocons=str(PHYLOCONS).lower()
    shell:
        "{params.juliacmd} {params.scriptpath:q} "
        "{output.genotypes} {params.outconsdir} {params.jlspath} "
        "{params.tree_groups:q} {params.known_genotypes:q} {params.catconsdir} "
        "{params.inconsdir} {params.blastdir} {params.minid} {params.phylocons}"

##################
# After checkpoint
##################
rule aln_groups:
    input: REFOUTDIR + "/groups/{segment}/{group}.fna"
    output: REFOUTDIR + "/tmp/aln_groups/{segment,[A-Z0-9]+}_{group}.aln.fna"
    log: "tmp/log/aln_groups/{segment}_{group}.aln.log"
    shell: "mafft {input:q} > {output} 2> {log}"

rule trim_groups:
    input: rules.aln_groups.output
    output: REFOUTDIR + "/groups/{segment}_{group}.aln.trim.fna"
    log: "tmp/log/aln_groups/{segment,[A-Z0-9]+}_{group}.trim.log"
    shell: "trimal -in {input} -out {output:q} -gt 0.9 -cons 60 2> {log}"

rule guide_tree:
    input: rules.trim_groups.output
    output: REFOUTDIR + "/group_trees/{segment}_{group}.treefile"
    log: "tmp/log/guide/{segment,[A-Z0-9]+}_{group}.log"
    threads: 2
    params: REFOUTDIR + "/group_trees/{segment}_{group}"
    shell: "iqtree -s {input:q} -pre {params:q} -T {threads} -m HKY+G2 --redo > {log}"

rule align_to_ref:
    input:
        ref=rules.trim_groups.output,
        con="tmp/catgroups/{segment}_{group}.fna"
    output: "tmp/merge/{segment,[A-Z0-9]+}_{group}.aln.fna"
    log: "tmp/log/merge/{segment}_{group}.log"
    shell: "mafft --add {input.con} --keeplength {input.ref:q} > {output} 2> {log}"

rule iqtree:
    input:
        aln=rules.align_to_ref.output,
        guide=rules.guide_tree.output
    output: "tmp/iqtree/{segment,[A-Z0-9]+}_{group}.treefile"
    log: "tmp/log/iqtree/{segment}_{group}.log"
    threads: 2
    params: "tmp/iqtree/{segment}_{group}"
    shell:
        "iqtree -s {input.aln} -pre {params} "
        "-g {input.guide:q} -T {threads} -m HKY+G2 --redo > {log}"

rule move_iqtree:
    input: rules.iqtree.output
    output: "trees/{segment,[A-Z0-9]+}/{group}.treefile"
    shell: "cp {input} {output}"

rule plot_tree:
    input: rules.move_iqtree.output
    output: "trees/{segment}/{group}.pdf"
    params: "tmp/catgroups/{segment}_{group}.fna"
    script: SNAKEDIR + "/scripts/plottree.py"
