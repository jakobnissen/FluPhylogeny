# Purpose:
# Given the internal tmp/internal.jls.gz file, create a CSV file with
# ngs run, sample name, good/bad, HA/NA type

using Consensus

function main(
    ngs_string::AbstractString,
    internal_path::AbstractString,
    genotype_path::AbstractString,
    outpath::AbstractString,
)
    parent = dirname(rstrip(rstrip(outpath, '/'), '\\'))
    (isempty(parent) || isdir(parent)) || error("Parent directory of output does not exist")
    ispath(outpath) && error("Output path exists: \"$outpath\"")

    isfile(internal_path) || error("No such path: \"$internal_path\"")
    isfile(genotype_path) || error("No such path: \"$genotype_path\"")

    internal = Consensus.load_internal(internal_path)
end
