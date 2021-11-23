# usage: julia gather_consensus.jl refdir tmp/cat consensus_dir

using FASTX: FASTA
using BioSequences: LongDNASeq
using Influenza: Segment, Sample, load_references, split_segment

"Represents a DNA sequence"
struct Seq
    name::String
    seq::LongDNASeq
end

function Seq(record::FASTA.Record)
    seq = FASTA.sequence(LongDNASeq, record)
    name = let
        h = FASTA.header(record)
        h === nothing ? error("Empty header in FASTA record") : h
    end
    Seq(name, seq)
end

function main(
    outdir::AbstractString,
    phylodir::AbstractString, # simply a dir with a list of segments
    consdir::AbstractString
)
    segments = map(readdir(phylodir)) do entry
        parse(Segment, entry)
    end

    consensus = load_consensus(consdir)

    # Remove consensus of irrelevant segments
    filter!(consensus) do (_, segment, _)
        in(segment, segments)
    end

    # Split by segment and ensure uniqueness of names within one segment
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

function load_consensus(
    dir::AbstractString
)::Vector{Tuple{Sample, Segment, Seq}}
    result = Vector{Tuple{Sample, Segment, Seq}}()
    record = FASTA.Record()
    for _sample in readdir(dir)
        sample = Sample(_sample)
        path = joinpath(dir, _sample, "primary.fna")
        isfile(path) && 
        open(FASTA.Reader, path) do reader
            while !eof(reader)
                read!(reader, record)
                header = let
                    h = FASTA.header(record)
                    h === nothing ? error("Record in $path has no header") : h
                end
                (name, segment) = split_segment(header)
                dnaseq = FASTA.sequence(LongDNASeq, record)
                seq = Seq(name, dnaseq)
                push!(result, (sample, segment, seq))
            end
        end
    end
    return result
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