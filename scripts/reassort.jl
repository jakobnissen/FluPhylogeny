module Reassort

using FASTX: FASTA
using InfluenzaCore: Segment
using BioSequences: LongDNASeq

const BlastRow = NamedTuple{
    (:query, :accession, :bitscore),
    Tuple{String, String, Float64}
}

ifilter(f) = x -> Iterators.filter(f, x)
imap(f) = x -> Iterators.map(f, x)

for T in (:Sample, :SubType)
    @eval begin
        struct $T
            name::String
        end
        Base.hash(x::$T, h::UInt8) = hash(x.name, h)
        Base.:(==)(x::$T, y::$T) = x.name == y.name
    end
end

function main(
    consensus_dir::AbstractString,
    blastdb::AbstractString
)
    check_blastn()
    fnaname, blastname = tempname(), tempname()
    namemap = copy_consensus(consensus_dir, fnaname)
    run(`blastn -db $blastdb -query $fnaname -outfmt '6 qacc sacc bitscore' -out $blastname -word_size 10`)
    blast_result = open(parse_blast, blastname)
    validate_blast!(blast_result, namemap)
end

function check_blastn()
    run(pipeline(`which blastn`, stdout=devnull, stderr=devnull))
end

function copy_consensus(
    consdir::AbstractString,
    outfile::AbstractString
)::Dict{String, Tuple{Sample, Segment}}
    namemap = Dict{String, Tuple{Sample, Segment}}()
    open(FASTA.Writer, outfile) do writer
        for subdir in readdir(consdir, join=true)
            sample = Sample(basename(subdir))
            open(FASTA.Reader, joinpath(subdir, "curated.fna")) do reader
                for record in reader
                    header = FASTA.header(record)
                    segment = parse(Segment, last(rsplit(header, '_', limit=2)))
                    haskey(namemap, header) && error("Sequence $header is not unique")
                    namemap[header] = (sample, segment)
                    write(writer, record)
                end
            end
        end
    end
    return namemap
end

function parse_blast(
    io::IO,
)::Vector{BlastRow}
    v = eachline(io) |>
    imap(strip) |>
    ifilter(!isempty) |>
    imap(i -> split(i, '\t')) |>
    imap() do (_query, _accession, _bitscore)
        query = String(_query)
        accession = String(_accession)
        bitscore = parse(Float64, _bitscore)
        (; query, accession, bitscore)
    end |> collect
end

function best!(v::Vector{BlastRow})
    isempty(v) && return v
    sort!(v, by=i -> (i.query, i.bitscore), rev=true)
    query = first(v).query
    len = 1
    for i in 2:lastindex(v)
        blastrow = v[i]
        if blastrow.query != query
            len += 1
            v[len] = blastrow
            query = blastrow.query
        end
    end
    resize!(v, len)
end

@noinline bad_trailing(s) = error("Cannot parse as NAME_SEGMENT: \"" * s, '"')

function split_segment(s::Union{String, SubString{String}})
    p = findlast(isequal(UInt8('_')), codeunits(s))
    p === nothing && return bad_trailing(s)
    seg = tryparse(Segment, SubString(s, p+1:lastindex(s)))
    seg === nothing && return bad_trailing(s)
    return (SubString(s, 1, prevind(s, p)), seg)
end

function validate_blast!(
    v::Vector{BlastRow},
    namemap::Dict{String, Tuple{Sample, Segment}}
)::Dict{Sample, Dict{Segment, Union{Nothing, SubType}}}
    result = Dict{Sample, Dict{Segment, Union{Nothing, SubType}}}()
    for i in best!(v)
        _, qseg = split_segment(i.query)
        _subtype, aseg = split_segment(i.accession)
        qseg == aseg || continue
        subtype = qseg == aseg ? SubType(_subtype) : nothing
        sample, segment = namemap[i.query]
        get!(valtype(result), result, sample)[segment] = subtype
    end
    return result
end

# TODO: Change this to give an output that makes sense
function report(subtypes::Dict{Sample, Dict{Segment, Union{Nothing, SubType}}})
    # Good: All segments map to same subtype
    good = Set((k for (k,v) in subtypes if !in(nothing, values(v)) && length(Set(values(v))) == 1))

    # Reassorted: Segments map to different subtypes
    reassorted = Set((k for (k,v) in setdiff(keys(subtypes), good) if length(setdiff(values(v), (nothing,))) > 1))

    # Bad: Neither good nor reassorted, meaning all mapping segments map to one subtype,
    # but some segments do not map
    bad = Set(setdiff(keys(subtypes), good, reassorted))
    return good, bad, reassorted
end

end # module