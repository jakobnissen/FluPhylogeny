# usage: julia gather_consensus.jl refdir tmp/catcons consensus_dir

module GatherConsensus

using FASTX: FASTA
using BioSequences: LongDNASeq
using Influenza: Segment, Sample, load_references, split_segment

struct SimpleFASTA
    header::String
    seq::LongDNASeq
end

"Represents a DNA sequence"
struct Seq
    name::String
    sample::Sample
    segment::Segment
    # order 1 means primary segment; 2,3,4 etc are secondary.
    # name + order combination is unique within a run
    order::UInt8
    seq::LongDNASeq
end

function main(
    outdir::AbstractString,
    phylodir::AbstractString, # simply a dir with a list of segments
    consdir::AbstractString
)
    segments = map(filter!(i -> !startswith(i, '.'), readdir(phylodir))) do entry
        parse(Segment, entry)
    end |> Set

    # Remove consensus of irrelevant segments
    consensus = filter!(seq -> seq.segment in segments, load_consensus(consdir))

    # Split by segment and ensure uniqueness of name+order within one segment
    bysegment = Dict(s => Seq[] for s in segments)
    name_orders = Dict(s => Set{Tuple{String, UInt8}}() for s in segments)
    for seq in consensus
        if (seq.name, seq.order) ∈ name_orders[seq.segment]
            error("Name \"$(seq.name)\", order $(seq.order) segment $(seq.segment) is not unique")
        end
        push!(bysegment[seq.segment], seq)
        push!(name_orders[seq.segment], (seq.name, seq.order))
    end

    dump_consensus(outdir, bysegment)
end

function load_consensus(dir::AbstractString)::Vector{Seq}
    result = Seq[]
    for _sample in filter!(i -> !startswith(i, '.'), readdir(dir))
        sample = Sample(_sample)
        append!(result, load_primary(joinpath(dir, _sample, "primary.fna"), sample))
        secondary_path = joinpath(dir, _sample, "secondary.fna")
        if isfile(secondary_path)
            append!(result, load_secondary(joinpath(dir, _sample, "secondary.fna"), sample))
        end
    end
    return result
end

function load_fna(path::AbstractString)::Vector{SimpleFASTA}
    open(FASTA.Reader, path) do reader
        map(reader) do record
            header = let
                h = FASTA.header(record)
                h === nothing ? error("Error: No header in record in \"$path\"") : h
            end
            SimpleFASTA(header, FASTA.sequence(LongDNASeq, record))
        end
    end
end

function load_primary(path::AbstractString, sample::Sample)::Vector{Seq}
    map(load_fna(path)) do fasta
        name, segment = split_segment(fasta.header)
        Seq(name, sample, segment, 1, fasta.seq)
    end
end

function load_secondary(path::AbstractString, sample::Sample)::Vector{Seq}
    map(load_fna(path)) do fasta
        name, segstr, order = rsplit(fasta.header, '_', limit=3)
        Seq(name, sample, parse(Segment, segstr), parse(UInt8, order), fasta.seq)
    end
end

function dump_consensus(
    outdir::AbstractString,
    consensus::Dict{Segment, Vector{Seq}}
)
    for (segment, seqs) in consensus
        open(FASTA.Writer, joinpath(outdir, string(segment) * ".fna")) do writer
            for seq in seqs
                name_order = seq.name * '_' * string(seq.order)
                write(writer, FASTA.Record(name_order, seq.seq))
            end
        end
    end 
end

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) == 3
        main(ARGS...)
    else
        println("Usage: julia gather_consensus.jl refdir tmp/cat consensus_dir")
        exit(1)
    end
end

end # module