module Tools

using BlastParse: BlastParse
using BioSequences: LongDNASeq
using FASTX: FASTA
using Influenza: Segment, Sample, split_segment

parse_blast_io(::IO) = nothing # to please linter
eval(BlastParse.gen_blastparse_code(
    (:qacc, :sacc, :qcovhsp, :pident, :bitscore),
    :parse_blast_io
))

"Represents a DNA sequence"
struct Seq
    name::String
    seq::LongDNASeq
end

function Seq(record::FASTA.Record)
    seq = FASTA.sequence(LongDNASeq, record)
    name = let
        h = FASTA.header(record)
        h === nothing ? error("Empty header in FASTA record") : h
    end
    Seq(name, seq)
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
        path = joinpath(dir, _sample, "primary.fna")
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

end # module