# TODO: ADD TEMP DECLARATIONS TO INTERMEDIATE FILES
# TODO: Better error message when given bad input paths
# TODO: Root in PDF file??
# TODO: Better quality PDF file - img size, colors etc.

import os
import pathlib

SNAKEDIR = os.path.dirname(workflow.snakefile)
JULIA_COMMAND = f"julia --startup-file=no --project={SNAKEDIR}"

###############
# Parse config
###############
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
possible_hosts = set(os.listdir(os.path.join(TOP_REF_DIR, "phylo")))
if HOST not in possible_hosts:
    raise KeyError(f"Directory for host {HOST} not found in {os.path.join(TOP_REF_DIR, 'phylo')}")

REFDIR = os.path.join(TOP_REF_DIR, "phylo", HOST)
REFOUTDIR = os.path.join(TOP_REF_DIR, "refout", "phylo", HOST)
pathlib.Path(REFOUTDIR).mkdir(parents=True, exist_ok=True)

# Get consensus dir
if "consensus" not in config:
    raise KeyError("You must supply consensus directory: '--config consensus=/path/to/consensus'")

CONSENSUS_DIR = os.path.abspath(os.path.expanduser(config["consensus"]))
if not os.path.isdir(CONSENSUS_DIR):
    raise NotADirectoryError(CONSENSUS_DIR)

ALL_SEGMENTS = sorted(os.listdir(REFDIR))
ALL_SUBTYPES = {s: [] for s in ALL_SEGMENTS}
for segment in ALL_SEGMENTS:
    for file in sorted(os.listdir(os.path.join(REFDIR, segment))):
        ALL_SUBTYPES[segment].append(os.path.splitext(file)[0])

###########################
# This is the function that triggers the checkpoint. By accessing the blastout
# checkpoint, the DAG "short-circuits" all the middle rules until the file
# has been completed.
def all_inputs(wildcards):
    # This is the function that gets all the segment/type combinations into the DAG
    trigger = checkpoints.blastall.get()
    with open(f"tmp/flutypes.txt") as file:
        combos = list(map(lambda s: s.partition('\t'), filter(None, map(str.strip, file))))
        return [f"trees/{s}/{f}.pdf" for (s,_,f) in combos]

rule all:
    input: all_inputs
    output: "commit.txt"
    params: SNAKEDIR
    shell: "git -C {params} rev-parse --short HEAD > {output}"

###########################
# Before checkpoint: BLASTn
###########################
# We use BLASTn to determine which files should be created
rule cat_ref:
    input: lambda wc: [os.path.join(REFDIR, segment, i + ".fna") for i in ALL_SUBTYPES[wc.segment]]
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
    shell: "makeblastdb -in {input} -dbtype nucl"

rule gather_consensus:
    output: expand("tmp/cat/{segment}.fna", segment=ALL_SEGMENTS)
    params:
        juliacmd=JULIA_COMMAND,
        scriptpath=f"{SNAKEDIR}/scripts/gather_consensus.jl",
        refdir=REFDIR,
        consensus_dir=CONSENSUS_DIR
    shell: "{params.juliacmd} {params.scriptpath} tmp/cat {params.refdir} {params.consensus_dir}"

rule blastn:
    input:
        nin=REFOUTDIR + "/{segment}.fna" + ".nin",
        nhr=REFOUTDIR + "/{segment}.fna" + ".nhr",
        nsq=REFOUTDIR + "/{segment}.fna" + ".nsq",
        fna="tmp/cat/{segment}.fna"
    output: "tmp/blast/{segment}.blastout"
    params:
        basename=rules.cat_ref.output
    shell: "blastn -query {input.fna} -db {params.basename} -outfmt '6 qacc sacc qlen length pident bitscore' > {output}"

# TODO: Somehow warn if no hits were found for a segment
# This also creates tmp/cat/{segment}_{flutype}.fna for all flutypes
# but this cannot appear in the graph at this point, so it's
# not part of the rule
rule parse_blast:
    input:
        blast=rules.blastn.output,
        cons=rules.gather_consensus.output
    output: "tmp/flutypes/{segment}.txt"
    params:
        juliacmd=JULIA_COMMAND,
        scriptpath=f"{SNAKEDIR}/scripts/parse_blast.jl",
        segment=lambda wc: wc.segment,
        catdir="tmp/cat"
    shell:
        "{params.juliacmd} {params.scriptpath} "
        "{params.segment} {input.cons} {input.blast} {params.catdir} {output}"

checkpoint blastall:
    input: expand("tmp/flutypes/{segment}.txt", segment=ALL_SEGMENTS)
    output: "tmp/flutypes.txt"
    run:
        combos = set()
        for inpath in input:
            segment, _ = os.path.splitext(os.path.basename(inpath))
            with open(inpath) as infile:
                for line in filter(None, map(str.strip, infile)):
                    seqname, flutype = line.split('\t')
                    combos.add((segment, flutype))

        with open(output[0], "w") as outfile:
            for (segment, flutype) in sorted(combos):
                print(segment, flutype, sep='\t', file=outfile)

##################
# After checkpoint
##################
rule aln_ref:
    input: ancient(REFDIR + "/{segment}/{flutype}.fna")
    output: "tmp/refaln/{segment}_{flutype}.aln.fna"
    log: "tmp/log/refaln/{segment}_{flutype}.aln.log"
    shell: "mafft {input} > {output} 2> {log}"

rule trim_ref:
    input: ancient(rules.aln_ref.output)
    output: REFOUTDIR + "/{segment}_{flutype}.aln.trim.fna"
    log: "tmp/log/refaln/{segment}_{flutype}.trim.log"
    shell: "trimal -in {input} -out {output} -gt 0.9 -cons 60 2> {log}"

rule guide_tree:
    input: ancient(rules.trim_ref.output)
    output: "tmp/guide/{segment}_{flutype}.treefile"
    log: "tmp/log/guide/{segment}_{flutype}.log"
    threads: 2
    params: "tmp/guide/{segment}_{flutype}"
    shell: "iqtree -s {input} -pre {params} -T {threads} -m HKY+G2 --redo > {log}"

rule move_guide_tree:
    input: ancient(rules.guide_tree.output)
    output: REFOUTDIR + "/{segment}_{flutype}.treefile"
    shell: "cp {input} {output}"

rule align_to_ref:
    input:
        ref=rules.trim_ref.output,
        con="tmp/cat/{segment}_{flutype}.fna"
    output: "tmp/merge/{segment}_{flutype}.aln.fna"
    log: "tmp/log/merge/{segment}_{flutype}.log"
    shell: "mafft --add {input.con} --keeplength {input.ref} > {output} 2> {log}"

rule iqtree:
    input:
        aln=rules.align_to_ref.output,
        guide=REFOUTDIR + "/{segment}_{flutype}.treefile"
    output: "tmp/iqtree/{segment}_{flutype}.treefile"
    log: "tmp/log/iqtree/{segment}_{flutype}.log"
    threads: 2
    params: "tmp/iqtree/{segment}_{flutype}"
    shell:
        "iqtree -s {input.aln} -pre {params} "
        "-g {input.guide} -T {threads} -m HKY+G2 --redo > {log}"

rule move_iqtree:
    input: rules.iqtree.output
    output: "trees/{segment}/{flutype}.treefile"
    shell: "cp {input} {output}"

rule plot_tree:
    input: rules.move_iqtree.output
    output: "trees/{segment}/{flutype}.pdf"
    script: SNAKEDIR + "/scripts/plottree.py"