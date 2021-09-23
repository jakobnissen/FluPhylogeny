
include("tools.jl")
using .Tools
using FASTX: FASTA

function main(outpath::AbstractString, consdir::AbstractString)
    open(FASTA.Writer, outpath) do writer
        for (sample, segment, seq) in Tools.load_consensus(consdir)
            newname = name(sample) * '_' * string(segment)
            write(writer, FASTA.Record(newname, seq.seq))
        end
    end
end

if length(ARGS) != 2
    println("Usage: julia cat_consensus.jl out.fna consensus_dir")
    exit(1)
else
    main(ARGS...)
end