using BioSequences: BioSequences, LongDNASeq, DNAMer, canonical
using MinHash: MinHash, MinHashSketch, MinHasher
using FASTX: FASTA
using InfluenzaCore: Segment
using ErrorTypes: ErrorTypes, Option, some, none, unwrap, @unwrap_or, is_error

function main(
    outdir::AbstractString,
    ref_dir::AbstractString,
    consensus_dir::AbstractString,
)
    # Load reference and consensus
    refs = load_ref_seqs(ref_dir)
    segments = Set(map(first, refs))
    cons = load_consensus(consensus_dir, segments)
    assert_unique_names(refs, cons)

    flutypes = get_flutypes(cons, refs)
    write_cons_cat(outdir, flutypes)
    return nothing
end

const MINHASHER = MinHasher(500)
const K = 8

function MinHash.update!(mh::MinHasher, seq::LongDNASeq)
    MinHash.update!(mh, (canonical(mer) for mer in BioSequences.each(DNAMer{K}, seq)))
end

"Represents a sample of a virus with a consensus sequence to be analyzed"
struct Sample
    name::String
end

"A \"subtype\" of influenza, some monophyletic clade of a segment"
struct FluType
    name::String
end

for T in (Sample, FluType)
    @eval begin
        name(x::$T) = x.name
        Base.hash(x::$T, h::UInt) = hash(name(x), hash($T, h))
        Base.:(==)(x::$T, y::$T) = name(x) == name(y)
        Base.print(io::IO, x::$T) = print(io, name(x))
    end
end

"Represents a DNA sequence"
struct Seq
    name::String
    seq::LongDNASeq
    sketch::MinHashSketch
end

function Seq(record::FASTA.Record)
    seq = FASTA.sequence(LongDNASeq, record)
    name = let
        n = FASTA.header(record)
        n === nothing ? error("Empty FASTA header in file") : n
    end
    empty!(MINHASHER)
    MinHash.update!(MINHASHER, seq)
    sketch = MinHashSketch(MINHASHER)
    Seq(name, seq, sketch)
end

#######################################################
# Load refs seqs
#######################################################
function load_ref_seqs(
    refseqdir::AbstractString
)::Vector{Tuple{Segment, FluType, Vector{Seq}}}
    segments = map(readdir(refseqdir)) do subdir
        s = tryparse(Segment, subdir)
        s === nothing ? error("Could not parse directory as segment: \"$subdir\"") : s
    end
    result = Vector{Tuple{Segment, FluType, Vector{Seq}}}()
    for segment in segments
        for path in readdir(joinpath(refseqdir, string(segment)), join=true)
            filename = basename(path)
            endswith(filename, ".fna") || error("Filename should end with .fna: $filename")
            flutype = FluType(filename[1:prevind(filename, ncodeunits(filename)-3)])
            push!(result, (segment, flutype, open(read_fasta, path)))
        end
    end
    return result
end

function read_fasta(io::IO)::Vector{Seq}
    record = FASTA.Record()
    result = Seq[]
    reader = FASTA.Reader(io)
    while !eof(reader)
        read!(reader, record)
        push!(result, Seq(record))
    end
    close(reader)
    return result
end

#######################################################
# Load consensus seqs
#######################################################
function load_consensus(
    consensusdir::AbstractString,
    segments::Set{Segment}
)::Vector{Tuple{Segment, Sample, Seq}}
    result = Vector{Tuple{Segment, Sample, Seq}}()
    for sampledir in readdir(consensusdir)
        sample = Sample(sampledir)
        filename = joinpath(consensusdir, sampledir, "curated.fna")
        for seq in open(read_fasta, filename)
            segment = let
                s = tryparse(Segment, last(rsplit(seq.name, '_', limit=2)))
                if s === nothing
                    error(
                        "Cannot parse file \"$filename\" header ",
                        "\"$(seq.name)\" as NAME_SEGMENT"
                    )
                else
                    s
                end
            end
            segment in segments || continue
            push!(result, (segment, sample, seq))
        end
    end
    return result
end

function assert_unique_names(
    refs::Vector{Tuple{Segment, FluType, Vector{Seq}}},
    cons::Vector{Tuple{Segment, Sample, Seq}}
)
    # They need only be unique within one segment
    names = Dict{Segment, Set{String}}()
    for (segment, _, seqs) in refs, seq in seqs
        if seq.name ∈ get!(Set{String}, names, segment)
            error("Duplicate name in refs: \"$(seq.name)\"")
        end
        push!(names[segment], seq.name)
    end
    for (segment, _, seq) in cons
        if seq.name ∈ get!(Set{String}, names, segment)
            error("Consensus name already present: \"$(seq.name)\"")
        end
        push!(names[segment], seq.name)
    end
    return nothing
end

#######################################################
# Find best flutypes
#######################################################
function get_flutypes(
    consensus::Vector{Tuple{Segment, Sample, Seq}},
    references::Vector{Tuple{Segment, FluType, Vector{Seq}}},
)::Dict{Tuple{Segment, FluType}, Vector{Seq}}
    result = Dict{Tuple{Segment, FluType}, Vector{Seq}}()

    for (cons_segment, sample, cons_seq) in consensus
        best_type = none(FluType)
        best_overlaps = 0
        for (ref_segment, flutype, ref_seqs) in references
            cons_segment == ref_segment || continue
            for ref_seq in ref_seqs
                overlaps = MinHash.intersectionlength(cons_seq.sketch, ref_seq.sketch)
                if overlaps > best_overlaps
                    best_overlaps = overlaps
                    best_type = some(flutype)
                end
            end
        end
        if is_error(best_type)
            @warn "Sample $(name(sample)) segment $(string(segment)) does not look like any subtype"
        else
            key = (cons_segment, unwrap(best_type))
            push!(get!(Vector{Seq}, result, key), cons_seq)
        end
    end
    return result
end

function write_cons_cat(
    directory::AbstractString,
    bytype::Dict{Tuple{Segment, FluType}, Vector{Seq}},
)
    # Write FASTA files
    for ((segment, flutype), seqs) in bytype
        open(joinpath(directory, "$(segment)_$(name(flutype)).fna"), "w") do io
            for seq in seqs
                println(io, '>', seq.name)
                println(io, seq.seq)
            end
        end
    end
    return nothing
end

if length(ARGS) != 3
    error("Usage: julia minhash.jl outdir refdir consensusdir")
else
    main(ARGS...)
end