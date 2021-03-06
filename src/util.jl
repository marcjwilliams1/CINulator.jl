
function copynumberfrequency(cells::Array{cancercellCN, 1})
    #return frequency of copy number mutations, note that we check if
    # copy number can be 0 even if this is not the case

    N = length(cells)
    Nchr = cells[1].genome.N

    CNstates = hcat(map(x -> gettotalcn(x), cells)...)'
    maxCN = maximum(CNstates)
    frequencymatrix = zeros(maxCN + 1, 2*Nchr)

    for i in 1:2*Nchr
        frequencymatrix[:, i] = counts(CNstates[:, i], 0:maxCN)
    end

    return frequencymatrix ./ N
end

function copynumberfrequencyA(cells::Array{cancercellCN, 1})
    #return frequency of copy number mutations, note that we check if
    # copy number can be 0 even if this is not the case

    N = length(cells)
    Nchr = cells[1].genome.N

    CNstates = hcat(map(x -> getAallele(x), cells)...)'
    maxCN = maximum(CNstates)
    frequencymatrix = zeros(maxCN + 1, 2*Nchr)

    for i in 1:2*Nchr
        frequencymatrix[:, i] = counts(CNstates[:, i], 0:maxCN)
    end

    return frequencymatrix ./ N
end

function copynumberfrequencyB(cells::Array{cancercellCN, 1})
    #return frequency of copy number mutations, note that we check if
    # copy number can be 0 even if this is not the case

    N = length(cells)
    Nchr = cells[1].genome.N

    CNstates = hcat(map(x -> getBallele(x), cells)...)'
    maxCN = maximum(CNstates)
    frequencymatrix = zeros(maxCN + 1, 2*Nchr)

    for i in 1:2*Nchr
        frequencymatrix[:, i] = counts(CNstates[:, i], 0:maxCN)
    end

    return frequencymatrix ./ N
end

function sitefrequency(cells::Array{cancercellCN, 1})
    N = length(cells)
    Nchr = cells[1].genome.N

    CNstates = hcat(map(x -> gettotalcn(x), cells)...)'
    maxCN = maximum(CNstates)
    frequency = Int64[]
    for i in 1:Nchr
        append!(frequency,
            counts(filter(x -> x != 2, CNstates[:, i]), 0:maxCN))
    end

    return filter!(x -> x > 0, frequency)
end


function sitefrequency(cells::Array{cancercellCN, 1}, chr)
    N = length(cells)
    Nchr = cells[1].genome.N

    CNstates = hcat(map(x -> gettotalcn(x), cells)...)'
    CNstateschr = CNstates[:, chr]
    maxCN = maximum(CNstates)
    frequency = counts(filter(x -> x != 2, CNstates), 0:maxCN)

    return filter!(x -> x > 0, frequency)
end

function celldataframe_locus(cell)
    chr = 1:cell.genome.N
    chrarm = ["1p", "1q", "2p", "2q", "3p", "3q", "4p", "4q","5p", "5q", "6p", "6q", "7p", "7q", "8p", "8q","9p", "9q", "10p", "10q", "11p", "11q", "12p", "12q","13p", "13q", "14p", "14q", "15p", "15q", "16p", "16q", "17p", "17q", "18p", "18q", "19p", "19q", "20p", "20q", "21p", "21q", "22p", "22q"]
    CNstates = gettotalcn(cell)
    DF = DataFrame(locus = chrarm[1:2*cell.genome.N],
                cell_id = cell.id,
                state = CNstates,
                A = getAallele(cell),
                B = getBallele(cell),
                fitness = cell.b - cell.d)
    DF[!, :state_AS_phased] = map((A,B) -> "$(A)|$(B)", DF[!,:A], DF[!,:B])
    return DF
end


function celldataframe_withpos(cell;
    chrlengths = [249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560, 59373566])
    chr = 1:cell.chromosomes.N
    startseg = fill(1, cell.chromosomes.N)
    endseg = chrlengths[1:cell.chromosomes.N]
    CNstates = cell.chromosomes.CN
    DF = DataFrame(chr = map(x -> string(x), chr),
                start = startseg,
                endseg = endseg,
                cell_id = cell.id,
                state = CNstates,
                fitness = cell.b - cell.d)
    return DF
end

function mergecelldataframe_pos(cells,
    chrlengths = [249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560, 59373566])
    return vcat(map(x -> celldataframe(x, chrlengths = chrlengths), cells)...)
end

function mergecelldataframe_locus(cells)
    return vcat(map(x -> celldataframe_locus(x), cells)...)
end

function copynumbernoise(df; celldist = Beta(1))
    segments = unique(df[:endseg])
    for cell_id in unique(df[:cell_id])
        prob = rand(celldist)
        for i in segments
            if rand() > prob
                df[(df[:cell_id] .== cell_id) .& (df[:endseg] .== i), :state] = 2
            end
        end
    end

    return df
end

function gettotalcn(cell)
    p = map(x -> x.p.tot, cell.genome.CN)
    q = map(x -> x.q.tot, cell.genome.CN)
    return [p q]'[:]
end

function gettotalcn(genome::Genome)
    p = map(x -> x.p.tot, genome.CN)
    q = map(x -> x.q.tot, genome.CN)
    return [p q]'[:]
end

function getAScn(cell)
    p = map(x -> "$(x.p.A)|$(x.p.B)", cell.genome.CN)
    q = map(x -> "$(x.q.A)|$(x.q.B)", cell.genome.CN)
    return [p q]'[:]
end

function getAallele(cell)
    p = map(x -> x.p.A, cell.genome.CN)
    q = map(x -> x.q.A, cell.genome.CN)
    return [p q]'[:]
end

function getBallele(cell)
    p = map(x -> x.p.B, cell.genome.CN)
    q = map(x -> x.q.B, cell.genome.CN)
    return [p q]'[:]
end

function meanfitness(cells)
    mean(map(x -> x.b, cells) .- map(x -> x.d, cells))
end

function minfitness(cells)
    minimum(map(x -> x.b, cells) .- map(x -> x.d, cells))
end

function maxfitness(cells)
    maximum(map(x -> x.b, cells) .- map(x -> x.d, cells))
end

function mediangenotype(cells)
    convert(Array{Int64, 1}, map(x -> round(x), median(hcat(map(x -> gettotalcn(x), cells)...), dims = 2))[:])
end

function meangenotype(cells)
    mean(hcat(map(x -> gettotalcn(x), cells)...), dims = 2)
end

function meanploidy(cells)
    mean(mean(hcat(map(x -> gettotalcn(x), cells)...), dims = 2))
end

function samplecells(cells, pct::Float64)
    pct < 1.0 || error("pct must be a value > 1.0, currently it is $(pct).")
    Nsamples = convert(Int64, round(pct * length(cells)))
    idx = sample(1:length(cells), Nsamples, replace = false)
    return cells[idx]
end

function samplecells(cells, Nsamples::Int64)
    Nsamples < length(cells) || error("Nsamples must be an integer < $(length(cells)), currently it is $(Nsamples).")
    idx = sample(1:length(cells), Nsamples, replace = false)
    return cells[idx]
end

function mutationdataframe(cell)
    DF = DataFrame(
                cell_id = cell.id,
                labels = cell.labelvec
                )
end

function mergemutationdataframe(cells)
    return vcat(map(x -> mutationdataframe(x), cells)...)
end

function show(io::IO, simresult::SimResult)
    println("##################################")
    println("Max population size: $(simresult.N[end])")
    println("Max time: $(simresult.t[end])")
    println("Mean fitness = $(round(meanfitness(simresult.cells), digits = 3)), Max fitness = $(round(maxfitness(simresult.cells), digits = 3)), Min fitness = $(round(minfitness(simresult.cells), digits = 3))")
    println("Mean genotype:")
    println("$(meangenotype(simresult.cells))")
    println("N divisions: $(simresult.Ndivision)")
    println("N dead cells: $(simresult.Ndeadcells)")
    println("Proportion cells die: $(simresult.Ndeadcells / simresult.Ndivision)")
    #println(cells[1])
    println("##################################")
end
