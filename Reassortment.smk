import os

SNAKEDIR = os.path.dirname(workflow.snakefile)
JULIA_COMMAND = f"julia --startup-file=no --project={SNAKEDIR}"

###############
# Parse config
###############
# Get refdir
if "ref" not in config:
    raise KeyError("You must supply reference directory: '--config ref=/path/to/ref'")

REFDIR = os.path.abspath(os.path.expanduser(config["ref"]))
if not os.path.isdir(REFDIR):
    raise NotADirectoryError(REFDIR)

REFOUTDIR = os.path.join(REFDIR, "refout")

# Get consensus dir
if "consensus" not in config:
    raise KeyError("You must supply consensus directory: '--config consensus=/path/to/consensus'")

CONSENSUS_DIR = os.path.abspath(os.path.expanduser(config["consensus"]))
if not os.path.isdir(CONSENSUS_DIR):
    raise NotADirectoryError(CONSENSUS_DIR)

###############
# Rules
###############
rule all:
    input: "subtypes.txt"

rule makeblastdb:
    input: REFDIR + "/humansubtypes.fna"
    output:
        nin=REFOUTDIR + "/{segment}" + ".nin",
        nhr=REFOUTDIR + "/{segment}" + ".nhr",
        nsq=REFOUTDIR + "/{segment}" + ".nsq"
    params:
        outbase = REFOUTDIR + "/humansubtypes"
    shell: "makeblastdb -in {input} -out {params.outbase} -dbtype nucl"

rule cat_cons:
    output: "tmp/reassort/all.fna"
    params:
        juliacmd=JULIA_COMMAND,
        scriptpath=f"{SNAKEDIR}/scripts/cat_consensus.jl",
        consensus_dir=CONSENSUS_DIR
    shell: "{params.juliacmd} {params.scriptpath} tmp/reassort/all.fna {params.consensus_dir}"

rule blastn:
    input:
        nin=REFOUTDIR + "/humansubtypes.nin",
        fna="tmp/reassort/all.fna"
    output: "tmp/reassort/all.blastout"
    params: REFOUTDIR + "/humansubtypes"
    shell: "blastn -query {input.fna} -db {params} -outfmt '6 qacc sacc qcovhsp pident bitscore' > {output}"

rule reassort:
    input:
        fasta="tmp/reassort/all.fna",
        blast=rules.blastn.output
    output: "subtypes.txt"
    params:
        juliacmd=JULIA_COMMAND,
        scriptpath=f"{SNAKEDIR}/scripts/reassort.jl",
    shell: "{params.juliacmd} {params.scriptpath} {input.blast} {input.fasta} {output}"