# usage: julia gather_consensus.jl refdir tmp/cat consensus_dir

using FASTX: FASTA
using InfluenzaCore: Segment

include("tools.jl")
using .Tools

function main(
    outdir::AbstractString,
    phylodir::AbstractString, # simply a dir with a list of segments
    consdir::AbstractString
)
    segments = map(readdir(phylodir)) do entry
        parse(Segment, entry)
    end

    consensus = Tools.load_consensus(consdir)

    # Remove consensus of irrelevant segments
    filter!(consensus) do (sample, segment, seq)
        in(segment, segments)
    end

    # Split by segment andensure uniqueness of names within one segment
    bysegment = Dict(s => Seq[] for s in segments)
    names = Dict(s => Set{String}() for s in segments)
    for (_, segment, seq) in consensus
        if seq.name ∈ names[segment]
            error("Name \"$(seq.name)\", segment $segment is not unique")
        end
        push!(bysegment[segment], seq)
        push!(names[segment], seq.name)
    end

    dump_consensus(outdir, bysegment)
end

function dump_consensus(
    outdir::AbstractString,
    consensus::Dict{Segment, Vector{Seq}}
)
    for (segment, seqs) in consensus
        open(FASTA.Writer, joinpath(outdir, string(segment) * ".fna")) do writer
            for seq in seqs
                write(writer, FASTA.Record(seq.name, seq.seq))
            end
        end
    end 
end

if length(ARGS) == 3
    main(ARGS...)
else
    println("Usage: julia gather_consensus.jl refdir tmp/cat consensus_dir")
    exit(1)
end