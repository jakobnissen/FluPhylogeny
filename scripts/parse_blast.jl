# Purpose -- applied to one segment at a time:
# Write "subtypes" with seqname \t flutype for all sequences
# Write tmp/cat/{segment}_{flutype}.fna for all present flutypes

include("tools.jl")
using .Tools
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
            seq = Seq(record)
            if !is_valid_seqname(seq.name)
                @warn "Invalid sequence name: \"$(seq.name)\". Will be renamed by IQ-TREE."
            end
            push!(seqs, seq)
        end
        return seqs
    end
end

function get_best(blastin::AbstractString)::Dict{String, FluType}
    rows = open(Tools.parse_blast_io, blastin)
    filter_blast!(rows)
    Tools.keep_best!(rows)
    return Dict(row.qacc => last(Tools.split_flutype(row.sacc)) for row in rows)
end

function filter_blast!(rows::Vector{<:NamedTuple})
    # For now, just some haphazardly chosen filters: Minimum 80% identity
    # over at least 80% of the query
    filter!(rows) do row
        row.pident ≥ 0.8 &&
        row.qcovhsp ≥ 0.8
    end
end

function write_best(
    segment::AbstractString,
    cat_dir::AbstractString,
    consensus::Vector{Seq},
    best::Dict{String, FluType}
)
    by_flutype = Dict{FluType, Vector{Seq}}()
    for seq in consensus
        if haskey(best, seq.name)
            flutype = best[seq.name]
            push!(get!(valtype(by_flutype), by_flutype, flutype), seq)
        else
            @warn "Sequence \"$(seq.name)\" does not match any flutype!"
        end
    end
    for (flutype, seqs) in by_flutype
        open(FASTA.Writer, joinpath(cat_dir, "$(segment)_$(name(flutype)).fna")) do writer
            for seq in seqs
                write(writer, FASTA.Record(seq.name, seq.seq))
            end
        end
    end
end

@noinline bad_trailing(s) = error("Cannot parse as NAME_FLUTYPE: \"" * s, '"')
function split_flutype(s_::Union{String, SubString{String}})
    s = strip(s_)
    p = findlast(isequal(UInt8('_')), codeunits(s))
    p === nothing && return bad_trailing(s)
    p == lastindex(s) && return bad_trailing(s)
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