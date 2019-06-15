
function copynumberfrequency(cells::Array{cancercellCN, 1})
    #return frequency of copy number mutations, note that we check if
    # copy number can be 0 even if this is not the case

    N = length(cells)
    Nchr = cells[1].chromosomes.N

    CNstates = hcat(map(x -> x.chromosomes.CN, cells)...)'
    maxCN = maximum(CNstates)
    frequencymatrix = zeros(maxCN + 1, Nchr)

    for i in 1:Nchr
        frequencymatrix[:, i] = counts(CNstates[:, i], 0:maxCN)
    end

    return frequencymatrix ./ N
end

function sitefrequency(cells::Array{cancercellCN, 1})
    N = length(cells)
    Nchr = cells[1].chromosomes.N

    CNstates = hcat(map(x -> x.chromosomes.CN, cells)...)'
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
    Nchr = cells[1].chromosomes.N

    CNstates = hcat(map(x -> x.chromosomes.CN, cells)...)'
    CNstateschr = CNstates[:, chr]
    maxCN = maximum(CNstates)
    frequency = counts(filter(x -> x != 2, CNstates), 0:maxCN)

    return filter!(x -> x > 0, frequency)
end


function celldataframe(cell; id = randstring(10),
    chrlengths = [249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560, 59373566])
    chr = 1:cell.chromosomes.N
    startseg = fill(1, cell.chromosomes.N)
    endseg = chrlengths[1:cell.chromosomes.N]
    CNstates = cell.chromosomes.CN
    DF = DataFrame(chr = chr,
                start = startseg,
                endseg = endseg,
                cell_id = id,
                state = CNstates,
                fitness = cell.b - cell.d)
    return DF
end

function mergecelldataframe(cells,
    chrlengths = [249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560, 59373566])
    return vcat(map(x -> celldataframe(x, chrlengths = chrlengths), cells)...)
end
