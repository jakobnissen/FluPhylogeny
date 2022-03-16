# Purpose: To concatenate
# 1) All reference segtypes to pan-segment references to find best segtype per segment
# 2) All references in tree groups to create guide trees and reference alns
module CatRef

using BioSequences: LongDNASeq
using FASTX: FASTA
using Influenza: Segment, Clade
using Phylogeny
using ErrorTypes: @unwrap_or

function main(
    refdir::AbstractString,
    known_genotypes_path::AbstractString,
    refoutdir::AbstractString
)
    known_genotypes = load_known_genotypes(known_genotypes_path)
    tree_groups = open(joinpath(refdir, "tree_groups.txt")) do io
        load_tree_groups(io, known_genotypes)
    end

    segtypes = load_segtypes(joinpath(refdir, "segments"), known_genotypes)
    dump_cat(joinpath(refoutdir, "segments"), segtypes)
    dump_groups(joinpath(refoutdir, "groups"), segtypes, tree_groups)

    touch(joinpath(refoutdir, "phylo_ref_cat_done"))
end

function load_segtypes(
    refdir::AbstractString,
    known_genotypes::Vector{GenoType}
)::Dict{Segment, Dict{Clade, Vector{FASTA.Record}}}
    known_clades = Dict(s => Set{Clade}() for s in instances(Segment))
    result = Dict{Segment, Dict{Clade, Vector{FASTA.Record}}}()
    for g in known_genotypes, (i, m_clade) in enumerate(g.v)
        segment = Segment(i - 1)
        clade = @unwrap_or m_clade continue
        push!(known_clades[segment], clade)
        get!(valtype(result), result, segment)[clade] = FASTA.Record[]
    end

    record = FASTA.Record()
    for subdir in filter!(i -> !startswith(i, '.'), readdir(refdir))
        segment = parse(Segment, subdir)
        joindir = joinpath(refdir, subdir)
        for file in filter!(i -> !startswith(i, '.'), readdir(joindir))
            m = match(r"^([A-Za-z0-9]+)\.fna$", file)
            if m === nothing
                error(
                    "In dir \"$joindir\", file \"$file\" does not conform to " *
                    "regex ^([A-Za-z0-9]+)\\.fna\$"
                )
            end
            clade = Clade(m[1]::AbstractString)
            result[segment][clade] = open(FASTA.Reader, joinpath(joindir, file)) do reader
                records = FASTA.Record[]
                names = Set{String}()
                while !eof(reader)
                    read!(reader, record)
                    if FASTA.hasdescription(record)
                        error("In file $file, FASTA record \"$(FASTA.header(record))\" contains whitespace")
                    end
                    id = FASTA.identifier(record)::String
                    if in(id, names)
                        error("In file $file, FASTA record \"$(FASTA.identifier(record))\" is not unique")
                    end
                    push!(names, id)
                    push!(records, FASTA.Record(id * "_" * string(clade), FASTA.sequence(LongDNASeq, record)))
                end
                records
            end
        end
    end
    return result
end

"Just dump a concatenation of all segtypes within a segment"
function dump_cat(
    outdir::AbstractString,
    segtypes::Dict{Segment, Dict{Clade, Vector{FASTA.Record}}}
)::Nothing
    isdir(outdir) || mkdir(outdir)
    for (segment, d) in segtypes
        open(FASTA.Writer, joinpath(outdir, string(segment) * ".fna")) do writer
            for (_, records) in d, record in records
                write(writer, record)
            end
        end
    end
end

function dump_groups(
    outdir::AbstractString,
    segtypes::Dict{Segment, Dict{Clade, Vector{FASTA.Record}}},
    tree_groups::Dict{Tuple{Segment, Clade}, Vector{String}}
)::Nothing
    by_group = Dict{Tuple{Segment, String}, Vector{FASTA.Record}}()
    for ((segment, clade), groups) in tree_groups, group in groups
        append!(get!(valtype(by_group), by_group, (segment, group)), segtypes[segment][clade])
    end

    isdir(outdir) || mkdir(outdir)
    for segment in keys(segtypes)
        subdir = joinpath(outdir, string(segment))
        isdir(subdir) || mkdir(subdir)
    end

    for ((segment, group), records) in by_group
        subdir = joinpath(outdir, string(segment))
        open(FASTA.Writer, joinpath(subdir, group * ".fna")) do writer
            foreach(i -> write(writer, i), records)
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) != 3
        println("Usage: julia cat_ref.jl refdir genotypes_path refoutdir")
        exit(1)
    end
    main(ARGS...)
end

end # module

