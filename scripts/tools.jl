module Tools

using BlastParse: BlastParse
using BioSequences: LongDNASeq
using FASTX: FASTA
using InfluenzaCore: Segment

parse_blast_io(::IO) = nothing # to please linter
eval(BlastParse.gen_blastparse_code(
    (:qacc, :sacc, :qcovhsp, :pident, :bitscore),
    :parse_blast_io
))

for T in (:Sample, :FluType)
    @eval begin
        struct $T
            name::String
        end
        name(x::$T) = x.name
        Base.hash(x::$T, h::UInt) = hash(name(x), hash($T, h))
        Base.:(==)(x::$T, y::$T) = name(x) == name(y)
        Base.print(io::IO, x::$T) = print(io, name(x))
        Base.isless(x::$T, y::$T) = isless(name(x), name(y))
    end
end

"Represents a DNA sequence"
struct Seq
    name::String
    seq::LongDNASeq
end

function keep_best!(v::Vector{<:NamedTuple})
    isempty(v) && return v
    sort!(v, by=i -> (i.qacc, i.bitscore), rev=true)
    query = first(v).qacc
    len = 1
    for i in 2:lastindex(v)
        blastrow = v[i]
        if blastrow.qacc != query
            len += 1
            v[len] = blastrow
            query = blastrow.qacc
        end
    end
    resize!(v, len)
end

function load_consensus(
    dir::AbstractString
)::Vector{Tuple{Sample, Segment, Seq}}
    result = Vector{Tuple{Sample, Segment, Seq}}()
    record = FASTA.Record()
    for _sample in readdir(dir)
        sample = Sample(_sample)
        path = joinpath(dir, _sample, "curated.fna")
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

@noinline bad_trailing(s) = error("Cannot parse as NAME_SEGMENT: \"" * s, '"')

function split_segment(s::Union{String, SubString{String}})
    p = findlast(isequal(UInt8('_')), codeunits(s))
    p === nothing && return bad_trailing(s)
    seg = tryparse(Segment, SubString(s, p+1:lastindex(s)))
    seg === nothing && return bad_trailing(s)
    return (SubString(s, 1, prevind(s, p)), seg)
end

export Seq,
    Sample,
    FluType,
    name

end # module