
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
