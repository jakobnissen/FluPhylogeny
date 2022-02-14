# TODO: ADD TEMP DECLARATIONS TO INTERMEDIATE FILES
# TODO: Better error message when given bad input paths
# TODO: Root in PDF file??
# TODO: Better quality PDF file - img size, colors etc.

import os

SNAKEDIR = os.path.dirname(workflow.snakefile)
JULIA_COMMAND = f"JULIA_LOAD_PATH='{SNAKEDIR}' julia --startup-file=no"

###############
# Parse config
###############
# Get refdir
if "ref" not in config:
    raise KeyError("You must supply reference directory: '--config ref=/path/to/ref'")

TOP_REF_DIR = os.path.abspath(os.path.expanduser(config["ref"]))
if not os.path.isdir(TOP_REF_DIR):
    raise NotADirectoryError(TOP_REF_DIR)

TOP_REFOUT_DIR = os.path.join(TOP_REF_DIR, "refout")

# Get host
if "host" not in config:
    raise KeyError("You must supply host: '--config host=human'")
HOST = config["host"]

# Check refdir has a host subdir
possible_hosts = set(os.listdir(TOP_REF_DIR))
if HOST not in possible_hosts:
    raise KeyError(f"Directory for host {HOST} not found in {TOP_REF_DIR}")

MINIMUM_IDENTITY = 0.96 if HOST == 'human' else 0.8

REFDIR = os.path.join(TOP_REF_DIR, HOST)
REFOUTDIR = os.path.join(TOP_REFOUT_DIR, HOST)

# Get consensus dir
if "consensus" not in config:
    raise KeyError("You must supply consensus directory: '--config consensus=/path/to/consensus'")

CONSENSUS_DIR = os.path.abspath(os.path.expanduser(config["consensus"]))
if not os.path.isdir(CONSENSUS_DIR):
    raise NotADirectoryError(CONSENSUS_DIR)

ALL_SEGMENTS = sorted(os.listdir(os.path.join(REFDIR, "segments")))

ALL_SEGTYPES = {s: [] for s in ALL_SEGMENTS}
for segment in ALL_SEGMENTS:
    for file in sorted(os.listdir(os.path.join(REFDIR, "segments", segment))):
        ALL_SEGTYPES[segment].append(os.path.splitext(file)[0])

# Get tree segments
with open(os.path.join(REFDIR, "tree_segments.csv")) as file:
    try:
        TREE_SEGMENTS = next(file).strip().split(",")
    except StopIteration:
        TREE_SEGMENTS = []
    if set(TREE_SEGMENTS) - set(ALL_SEGMENTS):
        raise ValueError("Tree segments must be a subset of all segments")

# Check the validity of "genotypes.tsv": No clade which is not in ALL_SEGTYPES
with open(os.path.join(REFDIR, "genotypes.tsv")) as file:
    header = next(file).strip()
    if not header.startswith("genotype\t"):
        raise ValueError("Expected header in genotypes.tsv reference file")
    headerfields = header.split("\t")

    if set(headerfields[1:]) - {"genotype"} != set(ALL_SEGMENTS):
        raise KeyError("Segments in genotypes.tsv do not match reference segments")

    for line in filter(None, map(str.strip, file)):
        fields = line.split('\t')
        if len(fields) != len(ALL_SEGMENTS) + 1:
            raise ValueError("Not all rows in genotypes.tsv has same length")

        for (segment, clade) in zip(headerfields[1:], fields[1:]):
            if clade not in ALL_SEGTYPES[segment]:
                raise KeyError(f"Clade {clade}, segment {segment} in genotypes.tsv has no reference")

###########################
# This is the function that triggers the checkpoint. By accessing the blastout
# checkpoint, the DAG "short-circuits" all the middle rules until the file
# has been completed.
def all_inputs(wildcards):
    # This is the function that gets all the segment/type combinations into the DAG
    trigger = checkpoints.genotypes.get()
    files = ["genotypes.txt", "tmp/genotypes.tsv"]
    
    combos = [p[:-4].partition('_') for p in os.listdir("tmp/cattypes")]
    if HOST == 'human':
        files.extend([f"trees/{s}/{f}.pdf" for (s,_,f) in combos])

    return files

rule all:
    input: all_inputs
    output: "phylogeny_versions.txt"
    params: SNAKEDIR
    shell:
        "cp {params:q}/copy_readme.md README_PHYLOGENY.md && "
        "git -C {params:q} rev-parse --short HEAD > {output} && "
        "julia -v >> {output}"

###########################
# Before checkpoint: BLASTn
###########################
rule instantiate:
    output: touch(TOP_REFOUT_DIR +  "/phylo_instantiated")
    params: JULIA_COMMAND
    shell: "{params} -e 'using Pkg; Pkg.resolve(); Pkg.instantiate()'"

# Cat all refs for each segment together as one to determine the best hit for each segment
rule cat_ref:
    input: lambda wc: [os.path.join(REFDIR, "segments", wc.segment, i + ".fna") for i in ALL_SEGTYPES[wc.segment]]
    output: REFOUTDIR + "/{segment}.fna"
    run:
        with open(output[0], "w") as outfile:
            for filename in input:
                with open(filename) as infile:
                    basename, _ = os.path.splitext(os.path.basename(filename))
                    for line in map(str.strip, infile):
                        if not line:
                            continue

                        if line.startswith('>'):
                            print(line + '_' + basename, file=outfile)
                        else:
                            print(line, file=outfile)

rule makeblastdb:
    input: rules.cat_ref.output
    output:
        nin=REFOUTDIR + "/{segment}.fna" + ".nin",
        nhr=REFOUTDIR + "/{segment}.fna" + ".nhr",
        nsq=REFOUTDIR + "/{segment}.fna" + ".nsq"
    log: "tmp/log/makeblastdb/{segment}.log"
    shell: "makeblastdb -in {input:q} -dbtype nucl 2> {log}"

rule gather_cons:
    input: TOP_REFOUT_DIR +  "/phylo_instantiated"
    output: expand("tmp/catcons/{segment}.fna", segment=ALL_SEGMENTS)
    params:
        juliacmd=JULIA_COMMAND,
        scriptpath=f"{SNAKEDIR}/scripts/gather_consensus.jl",
        segment_dir=os.path.join(REFDIR, "segments"),
        consensus_dir=CONSENSUS_DIR
    shell: "{params.juliacmd} {params.scriptpath:q} tmp/catcons {params.segment_dir:q} {params.consensus_dir:q}"

rule blastn:
    input:
        nin=REFOUTDIR + "/{segment}.fna" + ".nin",
        nhr=REFOUTDIR + "/{segment}.fna" + ".nhr",
        nsq=REFOUTDIR + "/{segment}.fna" + ".nsq",
        fna="tmp/catcons/{segment}.fna"
    output: "tmp/blast/{segment}.blastout"
    params:
        basename=rules.cat_ref.output
    shell: "blastn -query {input.fna} -db {params.basename} -outfmt '6 qacc sacc qcovhsp pident bitscore' > {output}"

# Purpose: Determine the subtype for each segment in ALL_SEGMENTS.
# write the output to tmp/clades.txt
# write human-readable report to genotypes.txt
checkpoint genotypes:
    input:
        blast=expand("tmp/blast/{segment}.blastout", segment=ALL_SEGMENTS),
        cons=expand("tmp/catcons/{segment}.fna", segment=ALL_SEGMENTS)
    output:
        clades="tmp/genotypes.tsv",
        genotypes="genotypes.txt"
        # Also tmp/cat/{segment}_{clade}.fna for every seg/type found
        # but the pairs ar unknown at this point, so it's not part of the rule
    params:
        juliacmd=JULIA_COMMAND,
        scriptpath=f"{SNAKEDIR}/scripts/parse_blast.jl",
        simple_genotypes=lambda wc: "simple.txt" if HOST == "swine" else "nothing",
        catconsdir="tmp/catcons", # corresponds to input.cons
        outconsdir="tmp/cattypes", # dir of output .fna files
        blastdir="tmp/blast", # corresponds to input.blast
        tree_segments=",".join(TREE_SEGMENTS),
        known_genotypes=os.path.join(REFDIR, "genotypes.tsv"),
        minid=MINIMUM_IDENTITY
    shell:
        "{params.juliacmd} {params.scriptpath:q} "
        "{output.genotypes} {params.simple_genotypes} {output.clades} {params.outconsdir} "
        "{params.tree_segments} {params.known_genotypes:q} {params.catconsdir} "
        "{params.blastdir} {params.minid}"

##################
# After checkpoint
##################

rule aln_ref:
    input: REFDIR + "/segments/{segment}/{clade}.fna"
    output: "tmp/refaln/{segment}_{clade}.aln.fna"
    log: "tmp/log/refaln/{segment}_{clade}.aln.log"
    shell: "mafft {input:q} > {output} 2> {log}"

rule trim_ref:
    input: rules.aln_ref.output
    output: REFOUTDIR + "/segments/{segment}_{clade}.aln.trim.fna"
    log: "tmp/log/refaln/{segment}_{clade}.trim.log"
    shell: "trimal -in {input} -out {output:q} -gt 0.9 -cons 60 2> {log}"

rule guide_tree:
    input: rules.trim_ref.output
    output: "tmp/guide/{segment}_{clade}.treefile"
    log: "tmp/log/guide/{segment}_{clade}.log"
    threads: 2
    params: "tmp/guide/{segment}_{clade}"
    shell: "iqtree -s {input:q} -pre {params} -T {threads} -m HKY+G2 --redo > {log}"

rule move_guide_tree:
    input: rules.guide_tree.output
    output: REFOUTDIR + "/{segment}_{clade}.treefile"
    shell: "cp {input} {output:q}"

rule align_to_ref:
    input:
        ref=rules.trim_ref.output,
        con="tmp/cattypes/{segment}_{clade}.fna"
    output: "tmp/merge/{segment}_{clade}.aln.fna"
    log: "tmp/log/merge/{segment}_{clade}.log"
    shell: "mafft --add {input.con} --keeplength {input.ref:q} > {output} 2> {log}"

rule iqtree:
    input:
        aln=rules.align_to_ref.output,
        guide=rules.move_guide_tree.output
    output: "tmp/iqtree/{segment}_{clade}.treefile"
    log: "tmp/log/iqtree/{segment}_{clade}.log"
    threads: 2
    params: "tmp/iqtree/{segment}_{clade}"
    shell:
        "iqtree -s {input.aln} -pre {params} "
        "-g {input.guide:q} -T {threads} -m HKY+G2 --redo > {log}"

rule move_iqtree:
    input: rules.iqtree.output
    output: "trees/{segment}/{clade}.treefile"
    shell: "cp {input} {output}"

rule plot_tree:
    input: rules.move_iqtree.output
    output: "trees/{segment}/{clade}.pdf"
    params: "tmp/cattypes/{segment}_{clade}.fna"
    script: SNAKEDIR + "/scripts/plottree.py"
