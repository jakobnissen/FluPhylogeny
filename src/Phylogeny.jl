module Phylogeny

using Influenza: Influenza, Clade, Sample, Segment, Segments, split_clade
using FASTX: FASTA
using ErrorTypes: Option, none, some, is_error, unwrap, @unwrap_or
using BioSequences: LongDNASeq

const N_SEGMENTS = length(instances(Segment))
const SegmentTuple{T} = NTuple{N_SEGMENTS, T}

ifilter(f) = x -> Iterators.filter(f, x)
imap(f) = x -> Iterators.map(f, x)

struct GenoType
    name::String
    # none if the segment is not considered
    v::SegmentTuple{Option{Clade}}
end

# Assumes format is .tsv
function load_known_genotypes(path::AbstractString)::Vector{GenoType}
    lines = eachline(path) |> imap(strip) |> ifilter(!isempty) |> imap(split) |> collect
    segments = map(i -> parse(Segment, i), lines[1][2:end])
    result = GenoType[]
    for line in lines[2:end]
        v = fill(none(Clade), N_SEGMENTS)
        for (segment, f) in zip(segments, line[2:end])
            v[Integer(segment) + 1] = some(Clade(f))
        end
        push!(result, GenoType(first(line), SegmentTuple(v)))
    end
    return sort!(result; by=i -> i.name)
end

function load_tree_groups(
    io::IO,
    known_genotypes::Vector{GenoType},
    # Map from (segment, clade) to list of all tree groups this clade is in
)::Dict{Tuple{Segment, Clade}, Vector{String}}
    tree_group_names = Dict(s => Set{String}() for s in instances(Segment))
    result = Dict{Tuple{Segment, Clade}, Vector{String}}()
    known_clades = Dict(s => Set{Clade}() for s in instances(Segment))
    for genotype in known_genotypes
        for (i, m_clade) in enumerate(genotype.v)
            segment = Segment(i - 1)
            clade = @unwrap_or m_clade continue
            push!(known_clades[segment], clade)
        end
    end
    header = split(readline(io)) # skip header
    if header != ["segment", "name", "segtypes"]
        error("Malformed header in tree_groups.txt")
    end
    for line in map(strip, eachline(io))
        isempty(line) && continue
        fields = split(line)
        if length(fields) != 3
            error(
                "In tree_groups.txt, expected 3 whitespace-delimited fields, got line:\n\"$line\"",
            )
        end
        segment = parse(Segment, fields[1])
        name = String(fields[2])
        if in(tree_group_names[segment], name)
            error("Duplicate tree group name: \"$name\" in $segment")
        end
        push!(tree_group_names[segment], name)
        for segtype in split(fields[3], ',')
            clade = Clade(segtype)
            if !in(clade, known_clades[segment])
                error(
                    "In tree_groups.txt, in group \"$name\", segtype" *
                    "\"$(segtype)\" is not a known segtype in genotypes.txt for segment $segment",
                )
            end
            push!(get!(valtype(result), result, (segment, clade)), name)
        end
    end
    return result
end

export GenoType,
    N_SEGMENTS, SegmentTuple, load_known_genotypes, load_tree_groups, ifilter, imap

end # module
