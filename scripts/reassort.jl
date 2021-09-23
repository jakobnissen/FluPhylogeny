include("tools.jl")
using .Tools
using InfluenzaCore: Segment

function main(
    blastpath::AbstractString,
    fastapath::AbstractString,
    outpath::AbstractString
)
    allpairs = open(get_all, fastapath)
    hits = open(read_blast, blastpath)
    open(outpath, "w") do io
        report(io, allpairs, hits)
    end
end

function get_all(io::IO)::Vector{Tuple{Sample, Segment}}
    result = Vector{Tuple{Sample, Segment}}()
    for line in eachline(io)
        if isempty(line) || first(line) != '>'
            continue
        end
        samplename, segment = Tools.split_segment(strip(line)[2:end])
        push!(result, (Sample(samplename), segment))
    end
    return result
end

function read_blast(io::IO)::Vector{Tuple{Sample, Segment, FluType}}
    rows = Tools.parse_blast_io(io)
    filter!(rows) do row
        row.pident ≥ 0.8 &&
        row.qcovhsp ≥ 0.8
    end
    Tools.keep_best!(rows)
    map(rows) do row
        sample, qsegment = Tools.split_segment(row.qacc)
        flutype, ssegment = Tools.split_segment(row.sacc)
        if qsegment != ssegment
            error("Sample $sample segment $qsegment maps to other segment $(row.sacc)")
        end
        (Sample(sample), qsegment, FluType(flutype))
    end
end

function report(
    io::IO,
    allpairs::Vector{Tuple{Sample, Segment}},
    hits::Vector{Tuple{Sample, Segment, FluType}}
)
    # Good: All segments map to same subtype
    # Reassorted: Segments map to different subtypes
    # Bad: Neither good nor reassorted, meaning all mapping segments map to one subtype,
    # but some segments do not map
    missings = let
        d = Dict{Sample, Set{Segment}}()
        for (sample, segment) in allpairs
            push!(get!(valtype(d), d, sample), segment)
        end
        for (sample, segment, _) in hits
            delete!(d[sample], segment)
        end
        filter!(d) do (k, v)
            !isempty(v)
        end
    end

    bysample = Dict{Sample, Vector{Tuple{Segment, FluType}}}()
    for (sample, segment, flutype) in hits
        push!(get!(valtype(bysample), bysample, sample), (segment, flutype))
    end

    good = Set((k for (k,v) in bysample if
        length(Set(map(last, v))) == 1 &&
        k ∉ keys(missings)
    ))

    reassorted = Set((k for (k,v) in bysample if
        length(Set(map(last, v))) > 1
    ))

    bad = Set(setdiff(keys(bysample), good, reassorted))

    if !isempty(good)
        println(io, "Good:")
        for sample in good
            println(io, '\t', sample, '\t', bysample[sample][1][2])
        end
    end

    if !isempty(bad)
        println(io, "\nSome segments unassigned:")
        for sample in bad
            println(io, '\t', sample, '\t', join(missings[sample], ", ", " and "))
        end
    end

    if !isempty(reassorted)
        println(io, "\nREASSORTMENT:")
        for sample in reassorted
            println(io, '\t', sample)
            for (segment, flutype) in sort!(bysample[sample])
                println(io, "\t\t", segment, '\t', flutype)
            end
        end
    end
end

if length(ARGS) != 3
    println("Usage: julia reassort.jl blastpath fastapath outpath")
    exit(1)
else
    main(ARGS...)
end