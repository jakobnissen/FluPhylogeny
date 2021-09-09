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
REFDIR = os.path.join(TOP_REF_DIR, "phylo", HOST)
if not os.path.isdir(REFDIR):
    raise NotADirectoryError(REFDIR)


REFOUTDIR = os.path.join(TOP_REF_DIR, "refout", "phylo", HOST)
pathlib.Path(REFOUTDIR).mkdir(parents=True, exist_ok=True)

# Get consensus dir
if "consensus" not in config:
    raise KeyError("You must supply consensus directory: '--config consensus=/path/to/consensus'")

CONSENSUS_DIR = os.path.abspath(os.path.expanduser(config["consensus"]))
if not os.path.isdir(CONSENSUS_DIR):
    raise NotADirectoryError(CONSENSUS_DIR)

###############
# Minhash
###############
# This script does the following:
# Read in all references, and passed consensus
# For each consensus, find the best suited flutype
# Cat all the consensus together in file called `tmp/cat/{SEGMENT}_{FLUTYPE}.fna`
# Write files in subtypes/{SEGMENT}.txt with the flutype for each sample. 
print("Determining closest subtypes")
if not os.path.exists("tmp"):
    os.mkdir("tmp")

if not os.path.exists("tmp/cat"):
    os.mkdir("tmp/cat")

if not os.path.exists("subtypes"):
    os.mkdir("subtypes")

shell(
    f"{JULIA_COMMAND} {os.path.join(SNAKEDIR, 'scripts/minhash.jl')} "
    "tmp/cat subtypes {REFDIR} {CONSENSUS_DIR}"
)

# Parse in the output
def get_pairs():
    files = list(filter(lambda f: f.endswith(".fna"), os.listdir("tmp/cat")))
    files = [f[:-4] for f in files]
    return [tuple(f.split("_")) for f in files]

# PAIRS is a list of (segment, flutype)
PAIRS = get_pairs()

###############
# Rules
###############
rule all:
    input: lambda wc: [f"trees/{s}/{t}.pdf" for (s, t) in PAIRS]

rule align_ref:
    input: ancient(REFDIR + "/{segment}/{flutype}.fna")
    output: "tmp/refaln/{segment}_{flutype}.aln.fna"
    log: "log/refaln/{segment}_{flutype}.aln.log"
    shell: "mafft {input} > {output} 2> {log}"

rule trim_aln:
    input: ancient(rules.align_ref.output)
    output: REFOUTDIR + "/{segment}/{flutype}.aln.trim.fna"
    log: "log/refaln/{segment}_{flutype}.trim.log"
    shell: "trimal -in {input} -out {output} -gt 0.9 -cons 60 2> {log}"

rule guide_tree:
    input: ancient(rules.trim_aln.output)
    output: "tmp/guide/{segment}_{flutype}.treefile"
    log: "log/guide/{segment}_{flutype}.log"
    threads: 2
    params: "tmp/guide/{segment}_{flutype}"
    shell: "iqtree -s {input} -pre {params} -T {threads} -m HKY+G2 --redo > {log}"

rule move_guide_tree:
    input: ancient(rules.guide_tree.output)
    output: REFOUTDIR + "/{segment}/{flutype}.treefile"
    shell: "cp {input} {output}"

rule merge_ref_cons:
    input:
        ref=rules.trim_aln.output,
        con="tmp/cat/{segment}_{flutype}.fna"
    output: "tmp/merge/{segment}_{flutype}.aln.fna"
    log: "log/merge/{segment}_{flutype}.log"
    shell: "mafft --add {input.con} --keeplength {input.ref} > {output} 2> {log}"

rule iqtree:
    input:
        aln=rules.merge_ref_cons.output,
        guide=REFOUTDIR + "/{segment}/{flutype}.treefile"
    output: "tmp/iqtree/{segment}_{flutype}.treefile"
    log: "log/iqtree/{segment}_{flutype}.log"
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