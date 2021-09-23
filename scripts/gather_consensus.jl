# usage: julia gather_consensus.jl refdir tmp/cat consensus_dir

using FASTX: FASTA
using InfluenzaCore: Segment
using BioSequences: LongDNASeq

function main(
    outdir::AbstractString,
    phylodir::AbstractString, # simply a dir with a list of segments
    consdir::AbstractString
)
    segments = map(readdir(phylodir)) do entry
        parse(Segment, entry)
    end

    consensus = load_consensus(segments, consdir)
    dump_consensus(outdir, consensus)
end

function load_consensus(
    segments::Vector{Segment},
    consdir::AbstractString
)::Dict{Segment, Dict{String, LongDNASeq}}
    result = Dict{Segment, Dict{String, LongDNASeq}}()
    record = FASTA.Record()
    for subdir in readdir(consdir, join=true)
        open(FASTA.Reader, joinpath(subdir, "curated.fna")) do reader
            while !eof(reader)
                read!(reader, record)
                header = FASTA.header(record)
                (name, segment) = let
                    s = try_split_segment(header)
                    if s === nothing
                        error("Header \"$header\" cannot be parsed as HEADER_SEGMENT")
                    end
                    s
                end
                in(segment, segments) || continue
                d = get!(valtype(result), result, segment)
                if haskey(d, name)
                    error("Name \"$name\", segment $segment is not unique")
                end
                d[name] = FASTA.sequence(LongDNASeq, record)
            end
        end
    end
    return result
end

function dump_consensus(
    outdir::AbstractString,
    consensus::Dict{Segment, Dict{String, LongDNASeq}}
)
    for (segment, d) in consensus
        open(FASTA.Writer, joinpath(outdir, string(segment) * ".fna")) do writer
            for (name, seq) in d
                write(writer, FASTA.Record(name, seq))
            end
        end
    end 
end

function try_split_segment(
    s_::Union{String, SubString{String}}
)::Union{Nothing, Tuple{SubString{String}, Segment}}
    s = strip(s_)
    p = findlast(isequal(UInt8('_')), codeunits(s))
    p === nothing && return nothing
    seg = tryparse(Segment, SubString(s, p+1, lastindex(s)))
    seg === nothing && return nothing
    return (SubString(s, 1:prevind(s, p)), seg)
end

if length(ARGS) == 3
    main(ARGS...)
else
    println("Usage: julia gather_consensus.jl refdir tmp/cat consensus_dir")
    exit(1)
end