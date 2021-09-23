using BlastParse: BlastParse
using BioSequences: LongDNASeq
using FASTX: FASTA

function main(
    segment::AbstractString,
    cons_path::AbstractString,
    blastin::AbstractString,
    cat_dir::AbstractString,
    subtypes::AbstractString
)
    consensus = load_consensus(cons_path)
    best = get_best(blastin)
    write_best(segment, cat_dir, consensus, best)

    open(subtypes, "w") do io
        for (seqname, flutype) in sort(collect(best))
            println(io, seqname, '\t', name(flutype))
        end
    end

    return nothing
end

eval(BlastParse.gen_blastparse_code(
    BlastParse.DEFAULT_COLUMNS,
    :parse_blast_io
))

"A \"flutype\" of influenza, some monophyletic clade of a segment"
struct FluType
    name::String
end

name(x::FluType) = x.name
Base.hash(x::FluType, h::UInt) = hash(name(x), hash(FluType, h))
Base.:(==)(x::FluType, y::FluType) = name(x) == name(y)
Base.print(io::IO, x::FluType) = print(io, name(x))
Base.isless(x::FluType, y::FluType) = isless(name(x), name(y))

"Represents a DNA sequence"
struct Seq
    name::String
    seq::LongDNASeq
end

function Seq(record::FASTA.Record)
    seq = FASTA.sequence(LongDNASeq, record)
    name = let
        n = FASTA.header(record)
        n === nothing ? error("Empty FASTA header in file") : n
    end
    if !is_valid_seqname(name)
        @warn "Invalid sequence name: \"$name\". Will be renamed by IQ-TREE."
    end
    Seq(name, seq)
end

"""IQ-TREE will rename any sequences that does not conform to this criteria.
It's probably better to warn here instead of later."""
function is_valid_seqname(s::AbstractString)
    all(s) do char
        isletter(char) ||
        char in '0':'9' ||
        char in ('.', '-', '_')
    end
end

# Input: The gathered consensus at tmp/cat/{segment}.fna
function load_consensus(cons_path::AbstractString)::Vector{Seq}
    record = FASTA.Record()
    open(FASTA.Reader, cons_path) do reader
        seqs = Seq[]
        while !eof(reader)
            read!(reader, record)
            push!(seqs, Seq(record))
        end
        return seqs
    end
end

function get_best(blastin::AbstractString)::Dict{String, FluType}
    bitscore = 0.0
    qacc = ""
    result = Dict{String, FluType}()
    rows = open(parse_blast_io, blastin)
    sort!(rows, by=row -> (row.qacc, row.bitscore), rev=true)
    for row in rows
        if row.qacc != qacc || row.bitscore > bitscore
            (_, flutype) = let
                s = try_split_flutype(row.sacc)
                if s === nothing
                    error("Cannot parse \"$(row.sacc)\" as NAME_FLUTYPE")
                else
                    s
                end
            end
            result[row.qacc] = flutype
            bitscore = row.bitscore
            qacc = row.qacc
        end
    end
    return result
end

function write_best(
    segment::AbstractString,
    cat_dir::AbstractString,
    consensus::Vector{Seq},
    best::Dict{String, FluType}
)
    by_flutype = Dict{FluType, Vector{Seq}}()
    for seq in consensus
        key = best[seq.name]
        push!(get!(valtype(by_flutype), by_flutype, key), seq)
    end

    for (flutype, seqs) in by_flutype
        open(FASTA.Writer, joinpath(cat_dir, "$(segment)_$(name(flutype)).fna")) do writer
            for seq in seqs
                write(writer, FASTA.Record(seq.name, seq.seq))
            end
        end
    end
end

function try_split_flutype(s_::Union{String, SubString{String}})
    s = strip(s_)
    p = findlast(isequal(UInt8('_')), codeunits(s))
    p === nothing && return nothing
    p == lastindex(s) && return nothing
    name = SubString(s, 1:prevind(s, p))
    flutype = FluType(SubString(s, p+1:lastindex(s)))
    return (name, flutype)
end

if length(ARGS) != 5
    println("Usage: julia parse_blast.jl segment catfna blastout catdir flutypes")
    exit(1)
else
    main(ARGS...)
end