# Purpose:
# Write "tmp/genotypes.tsv" with sample \t seqname \t clade for all sequences
# Create genotypes.txt report
# Write tmp/cat/{segment}_{clade}.fna for all present segment/clade combos

using Influenza
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
    v::Vector{Union{Clade, Nothing}}

    function GenoType(name::AbstractString, v::Vector{Union{Clade, Nothing}})
        if length(v) != N_SEGMENTS
            error("Must be $N_SEGMENTS long")
        end
        new(convert(String, name), v)
    end
end

@enum MatchFailure::UInt8 NoSegment NoMatch

function Base.print(io::IO, m::MatchFailure)
    m === NoSegment ? print(io, "No segment") :
    m === NoMatch ? print(io, "No match") :
    error()
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

struct SampleGenoType
    sample::Sample
    v::Vector{Union{Match, MatchFailure}}

    function SampleGenoType(name::Sample, v::Vector{Union{Match, MatchFailure}})
        if length(v) != N_SEGMENTS
            error("Must be $N_SEGMENTS long")
        end
        new(name, v)
    end
end

function clades(g::GenoType)
    ((Segment(i-1), f::Clade) for (i,f) in enumerate(g.v) if f isa Clade)
end

function clades(g::SampleGenoType)
    ((Segment(i-1), f.clade::Clade) for (i,f) in enumerate(g.v) if f isa Match)
end

function matches(g::SampleGenoType)
    ((Segment(i-1), f) for (i,f) in enumerate(g.v) if f isa Match)
end

function hits(g::SampleGenoType)
    ((Segment(i-1), f) for (i,f) in enumerate(g.v) if f != NoSegment)
end

# Could src be an instance of dst?
function is_compatible(sample::SampleGenoType, genotype::GenoType)
    all(1:N_SEGMENTS) do i
        s, g = sample.v[i], genotype.v[i]
        s == NoSegment || isnothing(g) || (s isa Match && s.clade == g)
    end
end

# is src actually an instance of dst?
function is_match(sample::SampleGenoType, genotype::GenoType)
    all(1:N_SEGMENTS) do i
        s, g = sample.v[i], genotype.v[i]
        isnothing(g) || (s isa Match && s.clade == g)
    end
end

function missing_segments(sample::SampleGenoType, genotype::GenoType)
    [Segment(i-1) for i in 1:N_SEGMENTS if sample.v[i] == NoSegment && genotype.v[i] isa Clade]
end

function main(
    genotype_report_path::AbstractString, # output: genotypes.txt
    simple_genotype_path::Union{Nothing, AbstractString},
    # output: Simplified version of genotypes.txt with only HA and NA, if present.
    genotypes_out_path::AbstractString,  # output: tmp/genotypes.tsv
    cattypes_dir::AbstractString, # output: dir to put {segment}_{clade}.fna
    tree_segments_str::AbstractString, # comma-sep string of segments to
        # write {segment}_{clade}.fna for
    known_genotypes_path::AbstractString, # input: genotypes.tsv ref input
    cat_dir::AbstractString, # input: dir w. concatenated consensus seqs
    blast_dir::AbstractString, # input: dir w. BLAST results
    minid::AbstractFloat
)
    isdir(cattypes_dir) || mkpath(cattypes_dir)
    consensus = load_consensus(cat_dir)
    sample_genotypes = load_sample_genotypes(blast_dir, consensus, minid)
    known_genotypes = load_known_genotypes(known_genotypes_path)

    # Output: tmp/genotypes.tsv
    write_genotypes(genotypes_out_path, simple_genotype_path, sample_genotypes)

    # Output: genotypes.txt
    write_genotype_report(genotype_report_path, known_genotypes, sample_genotypes)

    # Output: tmp/cat/{segment}_{clade}.fna for segments in TREE_SEGMENTS
    tree_segments = Set(map(i -> parse(Segment, i), split(tree_segments_str, ',')))
    
    write_fasta_combos(cattypes_dir, tree_segments, consensus, sample_genotypes)
    return nothing
end

function load_consensus(
    cons_dir::AbstractString # dir of {segment}.fna
)::Vector{Tuple{Sample, SegmentTuple{Option{LongDNASeq}}}}
    record = FASTA.Record()
    intermediate = Dict{Sample, Vector{Option{LongDNASeq}}}()
    for file in readdir(cons_dir, join=true)
        segment = parse(Segment, basename(first(splitext(file))))
        open(FASTA.Reader, file) do reader
            while !eof(reader)
                read!(reader, record)
                sample = let
                    h = FASTA.header(record)
                    h === nothing ? error("No header in record in $file") : Sample(h)
                end
                if !haskey(intermediate, sample)
                    intermediate[sample] = fill(none(LongDNASeq), N_SEGMENTS)
                end
                seq = FASTA.sequence(LongDNASeq, record)
                intermediate[sample][Integer(segment) + 1] = some(seq)
            end
        end
    end
    return [(s, SegmentTuple(v)) for (s, v) in intermediate]
end

function load_sample_genotypes(
    blast_dir::AbstractString,
    consensus::Vector{Tuple{Sample, SegmentTuple{Option{LongDNASeq}}}},
    minid::AbstractFloat,
)::Vector{SampleGenoType}
    
    # Initialize each segment with NoSegment
    sample_genotype_dict = Dict(
        sample => 
        fill!(Vector{Union{MatchFailure, Match}}(undef, N_SEGMENTS), NoSegment)
        for (sample, _) in consensus
    )

    # Set all the ones without missing segments to NoMatch
    for (sample, seqtuple) in consensus
        for (i, mseq) in enumerate(seqtuple)
            is_error(mseq) || (sample_genotype_dict[sample][i] = NoMatch)
        end
    end

    # Now read blast and set the matched segments to clades
    for file in readdir(blast_dir, join=true)
        segment = parse(Segment, basename(first(splitext(file))))
        rows = open(parse_blast_io, file)
        filter_blast!(rows, minid)
        keep_best!(rows)
        for row in rows
            sample = Sample(row.qacc)
            clade = last(split_clade(row.sacc))
            id = row.pident
            identifier = first(Influenza.split_clade(row.sacc))
            match = Match(clade, identifier, id)
            sample_genotype_dict[sample][Integer(segment) + 1] = match
        end
    end
    v = [SampleGenoType(k, v) for (k, v) in sample_genotype_dict]
    return sort!(v, by=i -> i.sample)
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

function filter_blast!(rows::Vector{<:NamedTuple}, minid::AbstractFloat)
    # For now, just some haphazardly chosen filters: Minimum 80% identity
    # over at least 80% of the query
    filter!(rows) do row
        row.pident ≥ minid &&
        row.qcovhsp ≥ 0.8
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
        push!(result, GenoType(first(line), v))
    end
    return sort!(result, by=i -> i.name)
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

# Output: genotypes.txt
function write_genotype_report(
    genotype_report_path::AbstractString,
    known_genotypes::Vector{GenoType},
    sample_genotypes::Vector{SampleGenoType}
)
    if isempty(known_genotypes) || isempty(sample_genotypes)
        touch(genotype_report_path)
        return nothing
    end

    good = Vector{Tuple{SampleGenoType, GenoType}}()
    indeterminate = Vector{Tuple{SampleGenoType, Vector{GenoType}}}()
    new = Vector{SampleGenoType}()

    possible_genotypes = GenoType[]
    for sample_genotype in sample_genotypes
        empty!(possible_genotypes)
        for genotype in known_genotypes
            if is_compatible(sample_genotype, genotype)
                push!(possible_genotypes, genotype)
            end
        end

        # If only one possibility, and precise match
        if length(possible_genotypes) == 1
            push!(good, (sample_genotype, first(possible_genotypes)))
        elseif isempty(possible_genotypes)
            push!(new, sample_genotype)
        else
            push!(indeterminate, (sample_genotype, copy(possible_genotypes)))
        end
    end

    if !isempty(new)
        @warn "Possible new genotype detected, see genotypes.txt report"
    end

    # Now write report
    open(genotype_report_path, "w") do io
        if !isempty(good)
            println(io, "Known genotypes:")
            for (s, g) in good
                print(io, '\t', s.sample, '\t', g.name)
                ms = missing_segments(s, g)
                if !isempty(ms)
                    print(io, " (missing: ", join(ms, ','), ')')
                end
                println(io)
            end
            println(io)
        end
        if !isempty(indeterminate)
            println(io, "Indeterminate genotypes:")
            for (g, gs) in indeterminate
                println(io, '\t', g.sample)
                for (seg, match) in matches(g)
                    sid = @sprintf "%.2f %%" (100 * match.identity)
                    println(io, "\t\t", join([seg, match.clade, match.identifier, sid], '\t'))
                end
                println(io)
                for g_ in gs
                    println(io, "\t\tMatches ", g_.name)
                end
                println(io)
            end
            println(io)
        end
        if !isempty(new)
            println(io, "New genotypes:")
            for g in new
                println(io, '\t', g.sample)
                for (seg, hit) in hits(g)
                    if hit == NoMatch
                        println(io, "\t\t", seg, '\t', hit)
                    else
                        sid = @sprintf "%.2f %%" (100 * hit.identity)
                        println(io, "\t\t", join([seg, hit.clade, hit.identifier, sid], '\t'))
                    end
                end
                println(io)
            end
        end
    end
    return nothing
end

function write_fasta_combos(
    cattypes_dir::AbstractString,
    tree_segments::Set{Segment},
    consensus::Vector{Tuple{Sample, SegmentTuple{Option{LongDNASeq}}}},
    sample_genotypes::Vector{SampleGenoType}
)
    genotype_of_sample = Dict(g.sample => g for g in sample_genotypes)
    by_combo = Dict{Tuple{Segment, Clade}, Vector{FASTA.Record}}()
    for (sample, mseqs) in consensus
        sample_genotype = genotype_of_sample[sample]
        for (segment, clade) in clades(sample_genotype)
            segment ∈ tree_segments || continue
            seq = @unwrap_or mseqs[Integer(segment) + 1] continue
            record = FASTA.Record(nameof(sample), seq)
            push!(get!(valtype(by_combo), by_combo, (segment, clade)), record)
        end
    end
    for ((segment, clade), records) in by_combo
        path = joinpath(cattypes_dir, "$(segment)_$(clade).fna")
        open(FASTA.Writer, path) do writer
            foreach(rec -> write(writer, rec), records)
        end
    end
    return nothing
end 

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) != 9
        println(
            "Usage: julia parse_blast.jl genotype_report_path " * 
            "genotype_out_path simple_genotype_path cattypes_dir tree_segments_str " *
            " known_genotypes_path cat_dir blast_dir min_id"
        )
        exit(1)
    else
        minid = parse(Float64, ARGS[9])
        if minid < 0 || minid > 1.0
            error("Minimum ID must be in 0..1")
        end
        simple_genotype_path = strip(ARGS[3]) == "nothing" ? nothing : ARGS[3]
        main(ARGS[1:2]..., simple_genotype_path, ARGS[4:8]..., minid)
    end
end
