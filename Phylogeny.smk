"""
Before DAG
One Julia call to:
    1. Load all possible refs
    2. Figure out which segments / samples passed.
    3. MinHash to determine the best for each

    Save as: 
        cat.fna for each (segment, flutype) of passed consensus

Then:
    * Parse the relevant pairs e.g. (segment, flutype) from cat filenames
    * Align + trim + guide tree each pair-ref
    * mafft add cat.fna to aligned pair-ref
    * iqtree added alignment
"""

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

# Get consensus dir
if "consensus" not in config:
    raise KeyError("You must supply consebsus directory: '--config consensus=/path/to/consensus'")

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
print("Determining closest subtypes")
if not os.path.exists("tmp"):
    os.mkdir("tmp")

if not os.path.exists("tmp/cat"):
    os.mkdir("tmp/cat")

shell(f"{JULIA_COMMAND} scripts/minhash.jl tmp/cat {REFDIR} {CONSENSUS_DIR}")

# Parse in the output
def get_pairs():
    files = list(filter(lambda f: f.endswith(".fna"), os.listdir("tmp/cat")))
    files = [f[:-4] for f in files]
    return [tuple(f.split("_")) for f in files]

# PAIRS is a list of (segment, flutype)
PAIRS = get_pairs()
print(PAIRS)

###############
# Rules
###############
# Step one: Get all flutypes aligned
# Step two: Build guide trees
# Step three: mafft --add cons to refaln
# Step four: iqtree above
# Step five: Visualization stuff