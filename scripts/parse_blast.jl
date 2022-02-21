# Purpose:
# Write "tmp/genotypes.tsv" with sample \t seqname \t clade for all sequences
# Create genotypes.txt report
# Write tmp/cat/{segment}_{clade}.fna for all present segment/clade combos

module ParseBlast

using Influenza: Influenza, Clade, Sample, Segment, Segments, split_clade
using FASTX: FASTA
using ErrorTypes: Option, none, some, is_error, unwrap, @unwrap_or
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
    # none if the segment is not considered
    v::SegmentTuple{Option{Clade}}
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
    # Top level none means no segment present. The key of the dict is their
    # "order" as assigned by consensus pipeline, to disambiguate between
    # segcopies
    v::SegmentTuple{Option{Dict{UInt8, Union{MatchFailure, Match}}}}
end

function is_empty(s::SampleGenoType, relevant::SegmentTuple{Bool})
    all(1:N_SEGMENTS) do i
        is_error(s.v[i]) || !relevant[i]
    end
end

function has_nohits_segment(s::SampleGenoType)
    any(s.v) do i
        d = @unwrap_or i return false
        any(isequal(NoHits), values(d))
    end
end

struct SampleSeqs
    sample::Sample
    # either missing segments, or any number of (order, seq)
    seqs::SegmentTuple{Option{Vector{Tuple{UInt8, LongDNASeq}}}}
end

@enum GenotypeMatch::UInt8 NotMatching CouldMatch UniqueMatch
function compatibility(sample::SampleGenoType, genotype::GenoType)::GenotypeMatch
    # If all segments are a Match to the right segtype, it's a unique match
    unique_match = true
    for i in 1:N_SEGMENTS
        s, g = sample.v[i], genotype.v[i]
        # Ignore segment if not considered by the Genotype
        genotype_clade = @unwrap_or g continue

        # If missing from sample, it can't be a unique match, but could still match
        dict = @unwrap_or s begin
            unique_match = false
            continue
        end

        # If any matches hit or any match is ambiguous, it could match.
        any_could = false
        for (_, maybe_match) in dict
            if maybe_match isa Match
                if maybe_match.clade == genotype_clade
                    any_could = true
                else
                    unique_match = false
                end
            elseif maybe_match == NoHits
                unique_match = false
            elseif maybe_match == AmbiguousHits
                unique_match = false
                any_could = true
            end
        end
        any_could || return NotMatching
    end
    return unique_match ? UniqueMatch : CouldMatch
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

    tree_segments = Set(map(i -> parse(Segment, i), split(strip(tree_segments_str), ',')))
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
    result = intermediate |> imap() do (sample, segment_orders)
        tup = ntuple(N_SEGMENTS) do i
            v = segment_orders[i]
            isempty(v) ? none : some(v)
        end
        SampleSeqs(sample, tup)
    end |> collect
    return sort!(result, by=i -> i.sample)
end

function load_sample_genotypes(
    blast_dir::AbstractString,
    consensus::Vector{SampleSeqs},
    min_id::AbstractFloat,
)::Vector{SampleGenoType}
    # Initialize with each present segment being set to no hits, since any segments
    # not present in the BLAST does indeed have no hits
    dType = Dict{UInt8, Union{MatchFailure, Match}}

    # TODO: Error in JET.jl: Try removing the "some": JET doesn't catch it!
    intermediate::Dict{Sample, SegmentTuple{Option{dType}}} = (consensus |> imap() do sampleseq
        tup::SegmentTuple{Option{dType}} = map(sampleseq.seqs) do maybe_segvector
            segvector = @unwrap_or maybe_segvector return none(dType)
            some(dType(o => NoHits for (o, s) in segvector))
        end
        sampleseq.sample => tup
    end |> Dict)

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
            d = unwrap(intermediate[sample][Integer(segment) + 0x01])
            d[order] = assign_clade(rowvec, min_id)
        end
    end 

    # Convert to output type
    result = intermediate |> imap() do (sample, st)
        SampleGenoType(sample, st)
    end |> collect 
    return sort!(result, by=i->i.sample)
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
        v = fill(none(Clade), N_SEGMENTS)
        for (segment, f) in zip(segments, line[2:end])
            v[Integer(segment) + 1] = some(Clade(f))
        end
        push!(result, GenoType(first(line), SegmentTuple(v)))
    end
    return sort!(result, by=i -> i.name)
end

"Get a SegmentTuple of whether a segment is relevant for assigning genotype"
function get_relevance(known_genotypes::Vector{GenoType})::SegmentTuple{Bool}
    is_relevant = fill(false, N_SEGMENTS)
    for genotype in known_genotypes
        v = genotype.v
        for i in 1:N_SEGMENTS
            if !is_error(v[i])
                is_relevant[i] = true
            end
        end
    end
    result = SegmentTuple(is_relevant)
    if all(!, result)
        error("No segment is relevant in known genotypes. Check reference genotypes.")
    end
    return result
end

# Categorize samples into empty, known genotypes, indetermine, and new
function categorize_genotypes(
    sample_genotypes::Vector{SampleGenoType},
    known_genotypes::Vector{GenoType}
)::NTuple{5, Vector}
    res_empty = SampleGenoType[]
    res_precise = Tuple{SampleGenoType, GenoType}[]
    # Not perfectly matching, but only one possibility of the known genotype
    res_unique = Tuple{SampleGenoType, GenoType}[]
    # unknown includes superinfection that cannot be nailed down to one genotype
    res_indeterminate = Tuple{SampleGenoType, Vector{GenoType}}[]
    res_new = SampleGenoType[]

    compatible_genotypes = GenoType[]
    matching_genotypes = GenoType[]

    relevant_segments = get_relevance(known_genotypes)

    for sample_genotype in sample_genotypes
        # Empty: All segments are either missing, or irrelevant for all genotypes.
        if is_empty(sample_genotype, relevant_segments)
            push!(res_empty, sample_genotype)
            continue
        end

        # If any segment has no hits, it's a new genotype
        if has_nohits_segment(sample_genotype)
            push!(res_new, sample_genotype)
            continue
        end

        # Else, we check compatibility
        empty!(compatible_genotypes)
        empty!(matching_genotypes)

        for genotype in known_genotypes
            match = compatibility(sample_genotype, genotype)
            if match == CouldMatch
                push!(compatible_genotypes, genotype)
            elseif match == UniqueMatch
                push!(matching_genotypes, genotype)
            else
                @assert match == NotMatching
            end
        end

        # Make sure the set of genotypes is not unresolvable
        if length(matching_genotypes) > 1
            error(
                "Sample \"$(sample_genotype.sample)\" matches multiple genotypes: " * 
                "\"$(matching_genotypes[1].name)\", and \"$(matching_genotypes[2].name)\""
            )
        end

        if !isempty(matching_genotypes) && !isempty(compatible_genotypes)
            error(
                "Sample \"$(sample_genotype.sample)\" matches genotype " * 
                "\"$(matching_genotypes[1].name)\", but is also compatible with " *
                "\"$(compatible_genotypes[1].name)\". Make sure genotypes are disjoint."
            )
        end

        # If there are any matching genotypes now, there can be only one, and no
        # compatible ones
        if !isempty(matching_genotypes)
            push!(res_precise, (sample_genotype, only(matching_genotypes)))
            continue
        end

        # If not compatible with anything, it's certainly a new one
        if isempty(compatible_genotypes)
            push!(res_new, sample_genotype)
            continue
        end

        # If only one genotype is compatible, that's the definition of uniquely matching
        if length(compatible_genotypes) == 1
            push!(res_unique, (sample_genotype, only(compatible_genotypes)))
            continue
        end

        # Else it's unknown
        push!(res_indeterminate, (sample_genotype, copy(compatible_genotypes)))
    end
    return (res_empty, res_precise, res_unique, res_indeterminate, res_new)
end

# Output: genotypes.txt
function write_genotype_report(
    genotype_report_path::AbstractString,
    res_empty::Vector{SampleGenoType},
    res_precise::Vector{Tuple{SampleGenoType, GenoType}},
    res_unique::Vector{Tuple{SampleGenoType, GenoType}},
    res_indeterminate::Vector{Tuple{SampleGenoType, Vector{GenoType}}},
    res_new::Vector{SampleGenoType}
)
    if !isempty(res_new)
        @warn "Possible new genotype detected, see genotypes.txt report"
    end
    open(genotype_report_path, "w") do io
        # Print empty
        if !isempty(res_empty)
            println(io, "Empty samples:")
            for sample_genotype in res_empty
                println(io, '\t', sample_genotype.sample)
            end
            println(io)
        end

        # Print known
        if !isempty(res_precise)
            println(io, "Precise genotypes:")
            for (sample_genotype, genotype) in res_precise
                println(io, '\t', sample_genotype.sample, '\t', genotype.name)
            end
            println(io)
        end

        # Print probable
        if !isempty(res_unique)
            println(io, "Uniquely matching genotypes:")
            for (sample_genotype, genotype) in res_unique
                println(io, '\t', sample_genotype.sample, '\t', genotype.name)
                print_match_status(io, sample_genotype)
                println(io)
            end
            println(io)
        end

        # Print unknown
        if !isempty(res_indeterminate)
            println(io, "Indeterminate genotypes:")
            for (sample_genotype, genotypes) in res_indeterminate
                println(io, '\t', sample_genotype.sample)
                print_match_status(io, sample_genotype)
                println(io)
                for genotype in genotypes
                    println(io, "\t\tMatches\t", genotype.name)
                end
                println(io)
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
        dict = @unwrap_or segment_status begin
            println(io, "\tNo segment")
            continue
        end 
        if length(dict) == 1
            m = only(values(dict))
            println(io, '\t', string(m))
        else
            println(io)
            for (order, m) in sort!(collect(dict), by=first)
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
            segment ∈ tree_segments || continue
            dict = @unwrap_or genotype.v[i] continue
            v = @unwrap_or sampleseq.seqs[i] continue
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