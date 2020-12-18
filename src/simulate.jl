mutable struct Arm
    A::Int64
    B::Int64
    tot::Int64
    function Arm(stateA, stateB)
        return new(stateA, stateB, stateA + stateB)
    end
end

mutable struct Chromosome
    p::Arm
    q::Arm
end

mutable struct Genome
    CN::Array{Chromosome, 1}
    N::Int64
    function Genome(N, states = [])
        if isempty(states)
            CN = fill(Chromosome(Arm(1, 1), Arm(1, 1)), N)
        else
            CN = map(x -> Chromosome(Arm(x[1]...), Arm(x[2]...)), states)
        end
        return new(CN, N)
    end
end

mutable struct cancercellCN
    b::Float64
    d::Float64
    binitial::Float64
    dinitial::Float64
    fitness::Array{Float64, 1}
    timedrivers::Array{Float64, 1}
    genome::Genome
    labelvec::Array{Int64, 1}
    id::String
end

mutable struct Chrmutrate
    mutrates::Array{Float64, 1}
    psplit::Array{Float64, 1}
    function Chrmutrate(N; m = 0.1, mutrates = [], psplit = [])
        if isempty(mutrates)
            mutrates = rand(Gamma(m), Int64(N))
        end
        if isempty(psplit)
            psplit = fill(0.01, N)
        end
        length(mutrates) == N || error("mutrates is length $(length(mutrates)), but number of chromosomes, N = $N. These should be equal")
        return new(mutrates, psplit)
    end
end

createalpha(alpha::Float64, N) = fill(alpha, N)
createalpha(alpha::Array{Float64, 1}, N) = alpha

mutable struct Chrfitness
    optimum::Genome
    alpha::Array{Float64, 1}
    function Chrfitness(optimumgenome::Genome; alpha = 0.0)
        return new(optimumgenome, createalpha(alpha, 2*length(optimumgenome.CN)))
    end
end

mutable struct Chrfitnessadditive
    s::Array{Float64, 1}
    function Chrfitnessadditive(N; s = [])
        if isempty(s)
            s = fill(1.0, N)
        end
        return new(s)
    end
end

mutable struct SimResult
    cells::Array{cancercellCN, 1}
    t::Array{Float64, 1}
    N::Array{Int64, 1}
    fitness::Array{Float64, 1}
    ploidy::Array{Float64, 1}
end

function copygenome(genomeold)
    genomenew = Genome(genomeold.N)
    for i in 1:length(genomenew.CN)
        genomenew.CN[i] = Chromosome(Arm(copy(genomeold.CN[i].p.A), copy(genomeold.CN[i].p.B)), Arm(copy(genomeold.CN[i].q.A), copy(genomeold.CN[i].q.B)))
        #genomenew.CN[i].q = Arm(copy(genomeold.CN[i].q.A), copy(genomeold.CN[i].q.B))
    end
    return genomenew
end

function copycell(cancercellold::cancercellCN)
  newcancercell::cancercellCN = cancercellCN(
  copy(cancercellold.b),
  copy(cancercellold.d),
  copy(cancercellold.binitial),
  copy(cancercellold.dinitial),
  copy(cancercellold.fitness),
  copy(cancercellold.timedrivers),
  copygenome(cancercellold.genome),
  copy(cancercellold.labelvec),
  cancercellold.id)
end

function initializesim(b, d, Nchr; N0 = 1, states = Genome(Nchr))

    #initialize time to zero
    t = 0.0
    tvec = Float64[]
    push!(tvec,t)

    N = N0
    Nvec = Int64[]
    push!(Nvec,N)

    #Initialize array of cell type that stores mutations for each cell and their fitness type
    #fitness type of 1 is the host population, lowest fitness
    cells = cancercellCN[]
    for i in 1:N0
        push!(cells, cancercellCN(b, d, b, d, Float64[], [], copygenome(states), Int64[], "temp"))
    end

    return t, tvec, N, Nvec, cells, [meanfitness(cells)], [meanploidy(cells)]
end

function additivefitness()
    
    f = function(cancercell, chrfitness, b, d)
    end

    return f
end

function genomedistance(genome1, genome2, alpha)

    dist = Float64[]
    for i in 1:length(genome1.CN)
        push!(dist, genome1.CN[i].p.tot - genome2.CN[i].p.tot)
        push!(dist, genome1.CN[i].q.tot - genome2.CN[i].q.tot)
    end

    return sqrt(sum(alpha .* dist.^2))
end

function optimumfitness(;increasebirth = true)
    if increasebirth == true
        f = function (cancercell, chrfitness, b, d)
            dist = genomedistance(cancercell.genome, chrfitness.optimum, chrfitness.alpha)
            smax = cancercell.binitial - cancercell.dinitial

            s = smax / (dist + 1)

            b = s + cancercell.dinitial
            d = b - s
            return b, d
        end
        return f
    else
        f = function (cancercell, chrfitness, b, d)
            dist = genomdistance(cancercell.genome, chrfitness.optimum, chrfitness.alpha)
            smax = cancercell.binitial - cancercell.dinitial

            s = smax / (dist + 1)

            d = maximum([cancercell.binitial - s, cancercell.dinitial])
            #d = cancercell.binitial - s
            b = s + d
            #d = cancercell.dinitial * (1 + propd * s)

            # b = (propb * b)  / (sum(abs.(distancefromoptimum)) * chrfitness.alpha + 1)
            # d = (propd * d)  * (sum(abs.(distancefromoptimum)) * chrfitness.alpha + 1)
            return b, d
        end
        return f
    end
end

function missegregation(cancercell1, cancercell2, chr, split_arm)

    #choose which cell and which allele are gained and lost
    cell_gain = rand(1:2)
    allele = sample(0:1, 2, replace = false)

    if (cancercell1.genome.CN[chr].p.A == 0) | (cancercell1.genome.CN[chr].q.A == 0)
        allele = [0, 1] 
    elseif (cancercell1.genome.CN[chr].p.B == 0) | (cancercell1.genome.CN[chr].q.B == 0)
        allele = [1, 0]
    end

    if split_arm == 1
        pq = rand(["p", "q"]) #randomly choose which chromosome arm is gained or lost
        if pq == "p"
            if cell_gain == 1
                cancercell1.genome.CN[chr].p = Arm(cancercell1.genome.CN[chr].p.A + allele[1], cancercell1.genome.CN[chr].p.B + allele[2])
                cancercell2.genome.CN[chr].p = Arm(cancercell2.genome.CN[chr].p.A - allele[1], cancercell2.genome.CN[chr].p.B - allele[2])
            else
                cancercell2.genome.CN[chr].p = Arm(cancercell2.genome.CN[chr].p.A + allele[1], cancercell2.genome.CN[chr].p.B + allele[2])
                cancercell1.genome.CN[chr].p = Arm(cancercell1.genome.CN[chr].p.A - allele[1], cancercell1.genome.CN[chr].p.B - allele[2])
            end
        else 
            if cell_gain == 1
                cancercell1.genome.CN[chr].q = Arm(cancercell1.genome.CN[chr].q.A + allele[1], cancercell1.genome.CN[chr].q.B + allele[2])
                cancercell2.genome.CN[chr].q = Arm(cancercell2.genome.CN[chr].q.A - allele[1], cancercell2.genome.CN[chr].q.B - allele[2])
            else
                cancercell2.genome.CN[chr].q = Arm(cancercell2.genome.CN[chr].q.A + allele[1], cancercell2.genome.CN[chr].q.B + allele[2])
                cancercell1.genome.CN[chr].q = Arm(cancercell1.genome.CN[chr].q.A - allele[1], cancercell1.genome.CN[chr].q.B - allele[2])
            end
        end

    else 

        if cell_gain == 1
            cancercell1.genome.CN[chr].p = Arm(cancercell1.genome.CN[chr].p.A + allele[1], cancercell1.genome.CN[chr].p.B + allele[2])
            cancercell1.genome.CN[chr].q = Arm(cancercell1.genome.CN[chr].q.A + allele[1], cancercell1.genome.CN[chr].q.B + allele[2])
            cancercell2.genome.CN[chr].p = Arm(cancercell2.genome.CN[chr].p.A - allele[1], cancercell2.genome.CN[chr].p.B - allele[2])
            cancercell2.genome.CN[chr].q = Arm(cancercell2.genome.CN[chr].q.A - allele[1], cancercell2.genome.CN[chr].q.B - allele[2])
        else
            cancercell2.genome.CN[chr].p = Arm(cancercell2.genome.CN[chr].p.A + allele[1], cancercell2.genome.CN[chr].p.B + allele[2])
            cancercell2.genome.CN[chr].q = Arm(cancercell2.genome.CN[chr].q.A + allele[1], cancercell2.genome.CN[chr].q.B + allele[2])
            cancercell1.genome.CN[chr].p = Arm(cancercell1.genome.CN[chr].p.A - allele[1], cancercell1.genome.CN[chr].p.B - allele[2])
            cancercell1.genome.CN[chr].q = Arm(cancercell1.genome.CN[chr].q.A - allele[1], cancercell1.genome.CN[chr].q.B - allele[2])
        end
        
    end

    return cancercell1, cancercell2
end

function newmutations(cancercell1, 
    cancercell2,
    μ::Chrmutrate,
    s::Chrfitness,
    Rmax, t,
    fitnessfunc;
    minCN = 1,
    maxCN = 6,
    labelcells = false,
    labelid = 1
    )

    mutations = map(x -> rand(Poisson(x)), μ.mutrates) .> 0
    #println(mutations)
    split_arms = map(x -> rand(Poisson(x)), μ.psplit) .> 0

    #record time if mutations occur
    if sum(mutations) > 0
        push!(cancercell1.timedrivers, t)
        push!(cancercell2.timedrivers, t)
    end

    killcell1 = false
    killcell2 = false

    #change copy number state of chromosomes
    for chr in 1:cancercell1.genome.N

        if mutations[chr] == 0
            continue
        end

        cancercell1, cancercell2 = missegregation(cancercell1, cancercell2, chr, split_arms[chr])

        #CN cannot go below minCN, cell dies
        if cancercell1.genome.CN[chr].p.tot < minCN || cancercell1.genome.CN[chr].q.tot < minCN
            #cancercell.chromosomes.CN[i] = minCN
            killcell1 = true
        end
        if cancercell2.genome.CN[chr].p.tot < minCN || cancercell2.genome.CN[chr].q.tot < minCN
            #cancercell.chromosomes.CN[i] = minCN
            killcell2 = true
        end
        #CN cannot exceed maxCN
        #if cancercell.genome.CN[i] > maxCN
        #    cancercell.genome.CN[i] = maxCN
        #end
    end


    b = cancercell1.binitial
    d = cancercell1.dinitial
    b, d = fitnessfunc(cancercell1, s, b, d)
    cancercell1.b = b
    cancercell1.d = d

    b = cancercell1.binitial
    d = cancercell1.dinitial
    b, d = fitnessfunc(cancercell2, s, b, d)
    cancercell2.b = b
    cancercell2.d = d

    if cancercell1.b + cancercell1.d > Rmax
      Rmax = cancercell1.b + cancercell1.d
    end

    if cancercell2.b + cancercell2.d > Rmax
        Rmax = cancercell2.b + cancercell2.d
      end

    if labelcells == true
        cancercell1.labelvec = push!(cancercell1.labelvec, labelid)
        labelid += 1
        cancercell2.labelvec = push!(cancercell2.labelvec, labelid)
        labelid += 1
    end

    return cancercell1, cancercell2, Rmax, killcell1, killcell2, labelid
end

exptime() = - log(rand())
meantime() = 1

function getfitness(cells, chrfitness, b, d; fitnessfunc = multiplicativefitness)
    for i in 1:length(cells)
        newb, newd = fitnessfunc(cells[i], chrfitness, b, d)
        cells[i].b = newb
        cells[i].d = newd
    end
    return cells
end

function simulate(b::Float64, d::Float64, Nmax::Int64, Nchr::Int64;
    N0 = 1, μ = Chrmutrate(Nchr, m = 0.01), s = Chrfitness(Genome(Nchr)),
    timefunction::Function = exptime, fitnessfunc = optimumfitness(),
    maxCN = 6, minCN = 1, states::Genome = Genome(Nchr),
    verbose = true,
    timestop = false, tend = 0.0, record = false,
    labelcells = false,
    nsamplecells = nothing)

    #Rmax starts with b + d and changes once a fitter mutant is introduced, this ensures that
    # b and d have correct units

    #initialize arrays and parameters
    t, tvec, N, Nvec, cells, fitness, ploidy = initializesim(b, d, Nchr, N0 = N0, states = states)
    labelid = 1
    cells = getfitness(cells, s, b, d, fitnessfunc = fitnessfunc)
    brate = maximum(map(x -> x.b, cells))
    drate = maximum(map(x -> x.d, cells))
    Rmax = brate + drate
    if tend == 0.0
        tend = (log(Nmax / N0) / (brate - drate)) + 0.571 / (brate - drate)
    end
    endsimulation = false
    if verbose
        println("##################################")
        println("Max population size: $(Nmax)")
        println("Max time: $(tend)")
        println("Rmax: $Rmax")
        println("Birth rate = $brate, death rate = $drate")
        println("initial Rmax: $Rmax")
        println("Mean fitness = $(meanfitness(cells)), Max fitness = $(maxfitness(cells)), Min fitness = $(minfitness(cells))")
        println("Initial distance from optimum: $(gettotalcn(cells[1]) .- gettotalcn(s.optimum))")
        #println(cells[1])
        println("##################################")
    end

    while endsimulation == false
        #pick a random cell
        randcell = rand(1:N)
        r = rand(Uniform(0, Rmax))
	    Nt = N
        Rmaxt = Rmax

        #birth event if r<birthrate, access correct birthrate from cells array

        #nothing if r > b+d
        if ((cells[randcell].b + cells[randcell].d) <= r )
          Δt =  1/(Rmax * Nt) * timefunction()
          t = t + Δt
        end

        #death event if b<r<b+d
        if (cells[randcell].b <= r < (cells[randcell].b + cells[randcell].d))
            #population decreases by 1
            N = N - 1
            #frequency of cell type decreases
            #remove deleted cell
            deleteat!(cells,randcell)
            Δt =  1/(Rmax * Nt) * timefunction()
            t = t + Δt
            #every cell dies reinitialize simulation
            if (N == 0)
                t, tvec, N, Nvec, cells, fitness, ploidy = initializesim(b, d, Nchr, N0 = N0, states = states)
            end
            continue
        end

        #birth event if b < r
        if r < cells[randcell].b

            #population increases by one
            N = N + 1
            #copy cell and mutations for cell that reproduces
            push!(cells,copycell(cells[randcell]))
            #add new mutations to both new cells
            cells[randcell], cells[end], Rmax, killcell1, killcell2, labelid =
                                            newmutations(cells[randcell], cells[end], μ,
                                            s, Rmax,t, fitnessfunc,
                                            maxCN = maxCN, minCN = minCN,
                                            labelcells = labelcells,
                                            labelid = labelid)
            if killcell1 == true
                N = N - 1
                deleteat!(cells,randcell)
            end
            if killcell2 == true
                N = N - 1
                deleteat!(cells,length(cells))
            end
            Δt =  1/(Rmaxt * Nt) * timefunction()
            t = t + Δt
        end

        #every cell dies reinitialize simulation
        if (N == 0)
            t, tvec, N, Nvec, cells, fitness, ploidy = initializesim(b, d, Nchr, N0 = N0, states = states)
        end

        if record
            push!(fitness, meanfitness(cells))
            push!(ploidy, meanploidy(cells))
        end
        push!(Nvec, N)
        push!(tvec,t)

        if ((timestop == true) & (t > tend))
            endsimulation = true
        elseif (timestop == false) & (N == Nmax)
            endsimulation = true
        end
        #println(CINulator.copynumberfrequencyA(cells))
    end
    if verbose
        println()
        println("##################################")
        println("Current time = $t, current N = $N")
        println("Mean fitness = $(meanfitness(cells)), Max fitness = $(maxfitness(cells)), Min fitness = $(minfitness(cells))")
        println("Median genotype:")
        println("$(mediangenotype(cells))")
        println("Mean genotype:")
        println("$(meangenotype(cells))")
        println("Optimum genotype:")
        println("$(gettotalcn(s.optimum))")
        println("Difference in genotype:")
        println("$(gettotalcn(s.optimum) .-mediangenotype(cells))")
        println("Average distance from optimum")
        println("$(sum(abs.(gettotalcn(s.optimum) .-meangenotype(cells))))")
        println("##################################")
        println()
    end

    for i in 1:length(cells)
        cells[i].id = randstring(10)
    end

    if nsamplecells != nothing
        nsamplecells = min(nsamplecells, Nmax)
        cells = samplecells(cells, nsamplecells)
    end

    return SimResult(cells, tvec, Nvec, fitness, ploidy)
end



function simulate(cells::Array{cancercellCN, 1}, tvec, Nvec, Nmax;
                    μ = Chrmutrate(Nchr, m = 0.01), 
                    s = Chrfitness(Genome(Nchr)),
                    timefunction::Function = exptime, 
                    fitnessfunc = optimumfitness(),
                    maxCN = 6, 
                    minCN = 1, 
                    states = [],
                    verbose = true,
                    timestop = false, 
                    tend = 0.0,
                    record = false,
                    nsamplecells = nothing)

    endsimulation = false

    sampledcells = deepcopy(cells)

    #Rmax starts with b + d and changes once a fitter mutant is introduced, this ensures that
    # b and d have correct units

    t = tvec[end]
    Nvec[end] = length(cells)
    N = Nvec[end]
    fitness = []
    ploidy = []

    #initialize arrays and parameters
    brate = maximum(map(x -> x.b, cells))
    drate = maximum(map(x -> x.d, cells))
    bratem = mean(map(x -> x.b, cells))
    dratem = mean(map(x -> x.d, cells))
    Rmax = brate + drate
    if tend == 0.0
        tend = (log(Nmax / N) / (bratem - dratem)) + 0.571 / (bratem - dratem)
    else
        tend = tend
    end

    ttic = 0.0

    if verbose
        println("##################################")
        println("Current time = $t, current N = $N")
        println("Birth rate = $brate, death rate = $drate")
        println("initial Rmax: $Rmax")
        println("Mean fitness = $(meanfitness(cells)), Max fitness = $(maxfitness(cells)), Min fitness = $(minfitness(cells))")
        println("Median genotype:")
        println("$(mediangenotype(cells))")
        println("Mean genotype:")
        println("$(meangenotype(cells))")
        println("Optimum genotype:")
        println("$(gettotalcn(s.optimum))")
        println("Difference in genotype:")
        println("$(gettotalcn(s.optimum) .-mediangenotype(cells))")
        println("Average distance from optimum")
        println("$(sum(abs.(gettotalcn(s.optimum) .-meangenotype(cells))))")
        #println(cells[1])
        println("##################################")
    end

    while endsimulation == false
        #pick a random cell
        randcell = rand(1:N)
        r = rand(Uniform(0, Rmax))
	    Nt = N
        Rmaxt = Rmax
        #birth event if r<birthrate, access correct birthrate from cells array

        #nothing if r > b+d
        if ((cells[randcell].b + cells[randcell].d) <= r )
          Δt =  1/(Rmax * Nt) * timefunction()
          t = t + Δt
          ttic = ttic + Δt
        end

        #death event if b<r<b+d
        if (cells[randcell].b <= r < (cells[randcell].b + cells[randcell].d))
            #population decreases by 1
            N = N - 1
            #frequency of cell type decreases
            #remove deleted cell
            deleteat!(cells,randcell)
            Δt =  1/(Rmax * Nt) * timefunction()
            t = t + Δt
            ttic = ttic + Δt
            #every cell dies reinitialize simulation
            if (N == 0)
                if verbose
                    println("Population size is 0, restarting")
                end
                cells = deepcopy(sampledcells)
                N = length(cells)
            end
            continue
        end

        #birth event if b < r
        if r < cells[randcell].b

            #population increases by one
            N = N + 1
            #copy cell and mutations for cell that reproduces
            push!(cells,copycell(cells[randcell]))
            #add new mutations to both new cells
            cells[randcell], Rmax, killcell, labelid =
                                                newmutations(cells[randcell], μ,
                                                s, Rmax, t, fitnessfunc,
                                                maxCN = maxCN,
                                                minCN = minCN,
                                                labelcells = labelcells,
                                                labelid = labelid)
            if killcell == true
                N = N - 1
                deleteat!(cells,randcell)
            end
            cells[end], Rmax, killcell = newmutations(cells[end], μ, s, Rmax, t,
                                            fitnessfunc, maxCN = maxCN,
                                            minCN = minCN,
                                            labelcells = labelcells,
                                            labelid = labelid)
            if killcell == true
                N = N - 1
                deleteat!(cells,length(cells))
            end
            Δt =  1/(Rmaxt * Nt) * timefunction()
            t = t + Δt
            ttic = ttic + Δt
        end

        #every cell dies reinitialize simulation
        if (N == 0)
            if verbose
                println("Population size is 0, restarting")
            end
            cells = deepcopy(sampledcells)
            N = length(cells)
        end

        if record
            push!(fitness, meanfitness(cells))
            push!(ploidy, meanploidy(cells))
        end
        push!(Nvec, N)
        push!(tvec,t)

        if ((timestop == true) & (ttic > tend))
            endsimulation = true
        elseif (timestop == false) & (N == Nmax)
            endsimulation = true
        end
    end
    if verbose
        println()
        println("##################################")
        println("Current time = $t, current N = $N")
        println("Mean fitness = $(meanfitness(cells)), Max fitness = $(maxfitness(cells)), Min fitness = $(minfitness(cells))")
        println("Median genotype:")
        println("$(mediangenotype(cells))")
        println("Mean genotype:")
        println("$(map(x -> round(x, digits = 3), meangenotype(cells)))")
        println("Optimum genotype:")
        println("$(gettotalcn(s.optimum))")
        println("Difference in genotype:")
        println("$(gettotalcn(s.optimum) .-mediangenotype(cells))")
        println("Average distance from optimum")
        println("$(map(x -> round(x, digits = 3), sum(abs.(gettotalcn(s.optimum).-meangenotype(cells)))))")
        println("##################################")
        println()
    end

    for i in 1:length(cells)
        cells[i].id = randstring(10)
    end

    return SimResult(cells, tvec, Nvec, fitness, ploidy)
end

function simulate_timeseries(b::Float64, d::Float64, Nmax::Int64, Nchr::Int64;
    N0 = 1, μ = Chrmutrate(Nchr, m = 0.01),
    s = Chrfitness(Nchr),
    timefunction::Function = exptime,
    fitnessfunc = optimumfitness(),
    maxCN = 6, minCN = 1, states = [],
    verbose = true, Ntimepoints = 5, pct = 0.1,
    timestop = false, tend = log(Nmax) / (b - d),
    record = false,
    nsamplecells = nothing)
    if verbose
        println("##################################")
        println("Timepoint 1")
        println("##################################")
        println()
    end
    simresult = simulate(b, d, Nmax::Int64, Nchr;
        N0 = N0, μ = μ, s = s,
        timefunction = timefunction,
        fitnessfunc = fitnessfunc,
        maxCN = maxCN,
        minCN = minCN,
        states = states,
        verbose = verbose,
        timestop = timestop,
        tend = tend,
        record = record,
        nsamplecells = nothing)
    simresult_t = []
    push!(simresult_t, simresult)
    sampledcells = samplecells(simresult.cells, pct)

    for i in 1:Ntimepoints-1
        if verbose
            println("##################################")
            println("Timepoint $(i + 1)")
            println("Number of cells sampled: $(length(sampledcells))")
            println("##################################")
            println()
        end
        simresult2 = simulate(sampledcells,
            simresult.t, simresult.N, Nmax;
            μ = μ, s = s,
            timefunction = timefunction,
            fitnessfunc = fitnessfunc,
            maxCN = maxCN,
            minCN = minCN,
            states = states,
            verbose = verbose,
            timestop = timestop,
            tend = tend,
            record = record,
            nsamplecells = nothing)
        push!(simresult_t, simresult2)
        sampledcells = samplecells(simresult2.cells, pct)
    end

    return simresult_t
end
