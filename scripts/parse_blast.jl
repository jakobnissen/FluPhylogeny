# Purpose:
# Write "tmp/genotypes.tsv" with sample \t seqname \t clade for all sequences
# Create genotypes.txt report
# Write tmp/cat/{segment}_{clade}.fna for all present segment/clade combos

module ParseBlast

using Influenza: Influenza, Clade, Sample, Segment, Segments, split_clade
using FASTX: FASTA
using ErrorTypes: Option, none, some, is_error, @unwrap_or
using BioSequences: LongDNASeq
using Printf: @sprintf
using BlastParse: BlastParse

const N_SEGMENTS = length(instances(Segment))
const SegmentTuple{T} = NTuple{N_SEGMENTS, T}

parse_blast_io(::IO) = nothing # to please linter
eval(BlastParse.gen_blastparse_code(
    (:qacc, :sacc, :qcovhsp, :pident, :bitscore),
    :parse_blast_io
))

ifilter(f) = x -> Iterators.filter(f, x)
imap(f) = x -> Iterators.map(f, x)

# Random optimization thought:
# If this ever bottlenecks, then enumerate all clades for each segment in 5 bits
# Store the genotype in a single int.
# Then use bitwise operations to compare genotypes
struct GenoType
    name::String
    # nothing if the segment is not considered
    v::SegmentTuple{Union{Clade, Nothing}}
end

struct Match
    clade::Clade
    identifier::String
    identity::Float64

    function Match(clade::Clade, identifier::AbstractString, identity::Number)
        id = convert(Float64, identity)
        (id ≤ 1.0 && id ≥ 0.0) || error("Identity must be in unit range")
        new(clade, convert(String, identifier), id)
    end
end

function Base.print(io::IO, m::Match)
    id = @sprintf "%.2f %%" (100 * m.identity)
    print(io, m.clade, '\t', m.identifier, '\t', id)
end

@enum MatchFailure::UInt8 NoHits AmbiguousHits

function Base.print(io::IO, m::MatchFailure)
    s = if m == NoHits
        "No match"
    elseif m == AmbiguousHits
        "Ambiguous BLAST hits"
    else
        error()
    end
    print(io, s)
end

struct SampleGenoType
    sample::Sample
    # Top level nothing means no segment present. Bottom-level nothing means
    # no match for the given segment. In the dict, the segments are ordered by
    # their "order"
    v::SegmentTuple{Union{Nothing, Dict{UInt8, Union{MatchFailure, Match}}}}
end

Base.isempty(s::SampleGenoType) = all(isnothing, s.v)

# This means: Does the sample genotype have multiple copies of the same
# segment that do not map to the same clade
function has_divergent_segments(s::SampleGenoType)
    any(s.v) do i
        # If a segment is missing, it is not divergent
        i === nothing && return false
        # If it's in only one copy, it's not divergent
        length(i) == 1 && return false
        clade = nothing
        any(i) do maybe_match
            # If there are multiple segments and one is not a match,
            # it's divergent
            maybe_match isa Match || return true
            if isnothing(clade)
                clade = match.clade
                false
            else
                # Alternatively, if it's a match but to a different
                # clade, it's divergent
                match.clade !== clade
            end
        end
    end
end

struct SampleSeqs
    sample::Sample
    # either missing segments, or any number of (order, seq)
    seqs::SegmentTuple{Union{Nothing, Vector{Tuple{UInt8, LongDNASeq}}}}
end

# Could the sample possibly be an instance of the genotype?
function is_compatible(sample::SampleGenoType, genotype::GenoType)::Bool
    all(1:N_SEGMENTS) do i
        s, g = sample.v[i], genotype.v[i]
        isnothing(g) || # specified segment is not considered in genotype
        isnothing(s) || # segment is missing from sample
        all(s) do (order, maybe_match)
            maybe_match isa Match &&
            maybe_match.clade == g
        end
    end
end

function main(
    genotype_report_path::AbstractString, # output: genotypes.txt
    cattypes_dir::AbstractString, # output: dir to put {segment}_{clade}.fna
    tree_segments_str::AbstractString, # comma-sep string of segments to
        # write {segment}_{clade}.fna for
    known_genotypes_path::AbstractString, # input: genotypes.tsv ref input
    cat_dir::AbstractString, # input: dir w. concatenated consensus seqs
    cons_dir::AbstractString, # input: dir to read a list of sample names from
    blast_dir::AbstractString, # input: dir w. BLAST results
    minid::AbstractFloat
)
    isdir(cattypes_dir) || mkpath(cattypes_dir)
    consensus = load_consensus(cat_dir, cons_dir)
    sample_genotypes = load_sample_genotypes(blast_dir, consensus, minid)
    known_genotypes = load_known_genotypes(known_genotypes_path)
    categories = categorize_genotypes(sample_genotypes, known_genotypes)
    write_genotype_report(genotype_report_path, categories...)

    tree_segments = Set(map(i -> parse(Segment, i), split(tree_segments_str, ',')))
    write_tree_fnas(cattypes_dir, tree_segments, consensus, sample_genotypes)
    return nothing
end

function load_consensus(
    cons_dir::AbstractString, # dir of {segment}.fna
    samples_dir::AbstractString # any dir to read samples from
)::Vector{SampleSeqs}
    record = FASTA.Record()
    intermediate = Dict(
        Sample(s) => 
        [Tuple{UInt8, LongDNASeq}[] for i in 1:N_SEGMENTS]
        for s in readdir(samples_dir)
    )
    for file in readdir(cons_dir, join=true)
        segment = parse(Segment, basename(first(splitext(file))))
        open(FASTA.Reader, file) do reader
            while !eof(reader)
                read!(reader, record)
                sample, order = let
                    h = FASTA.header(record)
                    h === nothing ? error("No header in record in $file") : Sample(h)
                    s, o = rsplit(h, '_')
                    (Sample(s), parse(UInt8, o))
                end
                seq = FASTA.sequence(LongDNASeq, record)
                push!(intermediate[sample][Integer(segment) + 1], (order, seq))
            end
        end
    end
    intermediate |> imap() do (sample, segment_orders)
        tup = ntuple(N_SEGMENTS) do i
            v = segment_orders[i]
            isempty(v) ? nothing : v
        end
        SampleSeqs(sample, tup)
    end |> collect
end

function load_sample_genotypes(
    blast_dir::AbstractString,
    consensus::Vector{SampleSeqs},
    min_id::AbstractFloat,
)::Vector{SampleGenoType}
    # The type inside SampleGenoType
    sT = SegmentTuple{Union{Nothing, Dict{UInt8, Union{MatchFailure, Match}}}}

    # Initialize with each present segment being set to no hits, since any segments
    # not present in the BLAST does indeed have no hits
    intermediate::Dict{Sample, sT} = Dict(map(consensus) do sampleseq
        tup = map(sampleseq.seqs) do maybe_segvector
            isnothing(maybe_segvector) && return nothing
            Dict(o => NoHits for (o, s) in maybe_segvector)
        end
        sampleseq.sample => tup
    end)

    # Now read blast and set the matched segments to clades
    for file in readdir(blast_dir, join=true)
        segment = parse(Segment, basename(first(splitext(file))))
        rows = open(parse_blast_io, file)
        byquery = Dict{Tuple{Sample, UInt8}, typeof(rows)}()
        for row in rows
            s, o = rsplit(row.qacc, '_', limit=2)
            sample, order = (Sample(s), parse(UInt8, o))
            push!(get!(valtype(byquery), byquery, (sample, order)), row)
        end
        for ((sample, order), rowvec) in byquery
            d = intermediate[sample][Integer(segment) + 0x01]::Dict
            d[order] = assign_clade(rowvec, min_id)
        end
    end 

    # Convert to output type
    intermediate |> imap() do (sample, st)
        SampleGenoType(sample, st)
    end |> collect 
end

# Pass in a vector that only contains one query
function assign_clade(
    v::Vector{<:NamedTuple},
    min_id::AbstractFloat
)::Union{Match, MatchFailure}
    # If the match covers less than 80% of the query, it doesn't matter
    filter!(i -> i.qcovhsp ≥ 0.8, v)

    # It's a match if: Hit above minimum identity
    # and second-best clade hit is 3%-points and 50% further away 
    sort!(v, by=i -> i.bitscore, rev=true)
    best = first(v)
    best.pident < min_id && return NoHits
    identifier, clade = split_clade(best.sacc)
    match = Match(clade, identifier, best.pident)
    next_clade_pos = findfirst(v) do row
        _, otherclade = split_clade(row.sacc)
        otherclade != clade
    end
    next_clade_pos === nothing && return match
    second_id = v[next_clade_pos].pident
    return if (second_id + 0.03 > best.pident ||
        (1 - second_id) < 1.5 * (1 - best.pident)
    )
        AmbiguousHits
    else
        match
    end
end

# Assumes format is .tsv
function load_known_genotypes(path::AbstractString)::Vector{GenoType}
    lines = eachline(path) |>
        imap(strip) |>
        ifilter(!isempty) |>
        imap(x -> split(x, '\t')) |>
        collect
    segments = map(i -> parse(Segment, i), lines[1][2:end])
    result = GenoType[]
    for line in lines[2:end]
        v = Vector{Union{Nothing, Clade}}(nothing, N_SEGMENTS)
        for (segment, f) in zip(segments, line[2:end])
            v[Integer(segment) + 1] = Clade(f)
        end
        push!(result, GenoType(first(line), SegmentTuple(v)))
    end
    return sort!(result, by=i -> i.name)
end

# The simple report is to make it easier for people to decide what to answer
# when a sample has been sent to NGS. Here, we just look at HA and NA
#=
function write_simple_report(
    path::AbstractString,
    sample_genotypes::Vector{SampleGenoType}
)
    open(path, "w") do io
        for genotype in sample_genotypes, (i, segment_result) in enumerate(genotype.v)
            segment = Segment(i - 1)
            # Status OK if 1 HA and NA clade

            # 

    end
end

function write_genotypes(
    path::AbstractString,
    simple_genotype_path::Union{Nothing, AbstractString},
    sample_genotypes::Vector{SampleGenoType}
)
    open(path, "w") do io
        println(io, "sample\tsegment\tclade")
        for genotype in sample_genotypes
            for (segment, clade) in clades(genotype)
                println(io, genotype.sample, '\t', segment, '\t', clade)
            end
        end
    end
    if simple_genotype_path !== nothing
        open(simple_genotype_path, "w") do io
            println(io, "sample\tgenotype")
            for genotype in sample_genotypes
                d = Dict(clades(genotype))
                print(io, genotype.sample, '\t')
                ha = get(d, Segments.HA, nothing)
                na = get(d, Segments.NA, nothing)
                str = if ha === nothing || na === nothing
                    "Unknown"
                else
                    string(ha) * string(na)
                end
                println(io, str)
            end
        end
    end
end
=#

# Categorize samples into empty, known genotypes, indetermine, and new
function categorize_genotypes(
    sample_genotypes::Vector{SampleGenoType},
    known_genotypes::Vector{GenoType}
)::NTuple{4, Vector}
    res_empty = SampleGenoType[]
    res_known = Tuple{SampleGenoType, GenoType}[]
    # unknown includes superinfection that cannot be nailed down to one genotype
    res_unknown = Tuple{SampleGenoType, Vector{GenoType}}[]
    res_new = SampleGenoType[]

    possible_genotypes = GenoType[]
    for sample_genotype in sample_genotypes
        if isempty(sample_genotype)
            push!(res_empty, sample_genotype)
            continue
        end

        empty!(possible_genotypes)
        for genotype in known_genotypes
            is_compatible(sample_genotype, genotype) && push!(possible_genotypes, genotype)
        end

        divergent = has_divergent_segments(sample_genotype)
        if length(possible_genotypes) == 1
            push!(res_known, (sample_genotype, only(possible_genotypes)))
        elseif length(possible_genotypes) == 0 && !divergent
            push!(res_new, sample_genotype)
        else
            push!(res_unknown, (sample_genotype, copy(possible_genotypes)))
        end
    end
    return (res_empty, res_known, res_unknown, res_new)
end

# Output: genotypes.txt
function write_genotype_report(
    genotype_report_path::AbstractString,
    res_empty::Vector{SampleGenoType},
    res_known::Vector{Tuple{SampleGenoType, GenoType}},
    res_unknown::Vector{Tuple{SampleGenoType, Vector{GenoType}}},
    res_new::Vector{SampleGenoType}
)
    if !isempty(res_new)
        @warn "Possible new genotype detected, see genotypes.txt report"
    end
    open(genotype_report_path, "w") do io
        # Print empty
        if !isempty(res_empty)
            println(io, "Empty genotypes:")
            for sample_genotype in res_empty
                println(io, '\t', sample_genotype.sample)
            end
            println(io)
        end

        # Print known
        if !isempty(res_known)
            println(io, "Known genotypes:")
            for (sample_genotype, genotype) in res_known
                println(io, '\t', sample_genotype.sample, '\t', genotype.name)
            end
            println(io)
        end

        # Print unknown
        if !isempty(res_unknown)
            println(io, "Indeterminate genotypes:")
            for (sample_genotype, genotypes) in res_unknown
                println(io, '\t', sample_genotype.sample)
                print_match_status(io, sample_genotype)
                println(io)
                for genotype in genotypes
                    println(io, "\t\tMatches\t", genotype.name)
                end
            end
            println(io)
        end

        # Print new
        if !isempty(res_new)
            println(io, "New genotypes:")
            for sample_genotype in res_new
                println(io, '\t', sample_genotype.sample)
                print_match_status(io, sample_genotype)
                println(io)
            end
            println(io)
        end
    end
end

function print_match_status(io::IO, s::SampleGenoType)
    for (i, segment_status) in enumerate(s.v)
        segment = Segment(i - 1)
        print(io, "\t\t", segment)
        if isnothing(segment_status)
            println(io, "\tNo segment")
        elseif length(segment_status) == 1
            m = only(values(segment_status))
            println(io, '\t', string(m))
        else
            println(io)
            for (order, m) in sort!(collect(segment_status), by=first)
                println(io, "\t\t\t", order, '\t', m)
            end
        end
    end
end

function write_tree_fnas(
    cattypes_dir::AbstractString,
    tree_segments::Set{Segment},
    sampleseqs::Vector{SampleSeqs},
    sample_genotypes::Vector{SampleGenoType}
)
    genotype_of_sample = Dict(g.sample => g for g in sample_genotypes)
    by_segtype = Dict{Tuple{Segment, Clade}, Vector{FASTA.Record}}()
    for sampleseq in sampleseqs
        genotype = genotype_of_sample[sampleseq.sample]
        for i in 1:N_SEGMENTS
            segment = Segment(i - 1)
            dict = genotype.v[i]
            dict === nothing && continue
            v = sampleseq.seqs[i]
            v === nothing && continue
            for (order, seq) in v
                maybe_match = dict[order]
                maybe_match isa Match || continue
                header = string(sampleseq.sample) * '_' * string(order)
                record = FASTA.Record(header, seq)
                key = (segment, maybe_match.clade)
                push!(get!(valtype(by_segtype), by_segtype, key), record)
            end
        end
    end

    for ((segment, clade), records) in by_segtype
        path = joinpath(cattypes_dir, "$(segment)_$(clade).fna")
        open(FASTA.Writer, path) do writer
            foreach(rec -> write(writer, rec), records)
        end
    end
    return nothing
end 

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) != 8
        println(
            "Usage: julia parse_blast.jl genotype_report_path " * 
            "cattypes_dir tree_segments_str " *
            "known_genotypes_path cat_dir cons_dir blast_dir min_id"
        )
        exit(1)
    else
        minid = parse(Float64, ARGS[8])
        if !isfinite(minid) || minid < 0 || minid > 1
            error("Minimum ID must be in 0..1")
        end
        main(ARGS[1:7]..., minid)
    end
end

end # module