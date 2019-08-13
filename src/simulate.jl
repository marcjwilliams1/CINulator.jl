mutable struct Chromosomes
    CN::Array{Int64, 1}
    N::Int64
    function Chromosomes(N, states = [])
        if isempty(states)
            CN = fill(2, N)
        else
            CN = states
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
    chromosomes::Chromosomes
end

mutable struct Chrmutrate
    loss::Array{Float64, 1}
    gain::Array{Float64, 1}
    function Chrmutrate(N; m = 0.1, mutratesgain = [],
        mutratesloss = [])
        if isempty(mutratesloss)
            mutratesloss = rand(Gamma(m), Int64(N))
        end
        if isempty(mutratesgain)
            mutratesgain = rand(Gamma(m), Int64(N))
        end
        length(mutratesgain) == N || error("mutratesgain is length $(length(mutratesgain)), but number of chromosomes, N = $N. These should be equal")
        length(mutratesloss) == N || error("mutratesloss is length $(length(mutratesloss)), but number of chromosomes, N = $N. These should be equal")
        return new(mutratesloss, mutratesgain)
    end
end

createalpha(alpha::Float64, N) = fill(alpha, N)
createalpha(alpha::Array{Float64, 1}, N) = alpha

mutable struct Chrfitness
    optimum::Array{Int64, 1}
    alpha::Array{Float64, 1}
    function Chrfitness(N; optimum = [], alpha = 0.0)
        if isempty(optimum)
            optimum = rand(1:6, N)
        end
        length(optimum) == N || error("optimum is length $(length(optimum)), but number of chromosomes, N = $N. These should be equal")
        return new(optimum, createalpha(alpha, N))
    end
end

function copychr(chromosomesnew, chromosomesold)
    chromosomesnew.CN = copy(chromosomesold.CN)
    return chromosomesnew
end

function copycell(cancercellold::cancercellCN)
  newcancercell::cancercellCN = cancercellCN(
  copy(cancercellold.b),
  copy(cancercellold.d),
  copy(cancercellold.binitial),
  copy(cancercellold.dinitial),
  copy(cancercellold.fitness),
  copy(cancercellold.timedrivers),
  copychr(Chromosomes(cancercellold.chromosomes.N),
  cancercellold.chromosomes))
end

function initializesim(b, d, Nchr; N0 = 1, states = [])

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
        push!(cells, cancercellCN(b, d, b, d, Float64[], [], Chromosomes(Nchr, states)))
    end

    return t, tvec, N, Nvec, cells
end

function optimumfitness(;increasebirth = true)
    if increasebirth == true
        f = function (cancercell, chrfitness, b, d)
            distancefromoptimum = cancercell.chromosomes.CN .- chrfitness.optimum
            dist = sum(chrfitness.alpha .* abs.(distancefromoptimum))
            smax = cancercell.binitial - cancercell.dinitial

            s = smax / (dist + 1)

            b = s + cancercell.dinitial
            d = b - s
            #d = cancercell.dinitial * (1 + propd * s)

            # b = (propb * b)  / (sum(abs.(distancefromoptimum)) * chrfitness.alpha + 1)
            # d = (propd * d)  * (sum(abs.(distancefromoptimum)) * chrfitness.alpha + 1)
            return b, d
        end
        return f
    else
        f = function (cancercell, chrfitness, b, d)
            distancefromoptimum = cancercell.chromosomes.CN .- chrfitness.optimum

            dist = sum(chrfitness.alpha .* abs.(distancefromoptimum))
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

function newmutations(cancercell::cancercellCN,
    μ::Chrmutrate,
    s::Chrfitness,
    Rmax, t,
    fitnessfunc;
    minCN = 1,
    maxCN = 6)

    #function to add new mutations to cells based on μ
    mutations_gain = map(x -> rand(Poisson(x)), μ.gain) .> 0
    mutations_loss = map(x -> rand(Poisson(x)), μ.loss) .> 0

    #record time if mutations occur
    if (sum(mutations_gain) + sum(mutations_loss)) > 0
        push!(cancercell.timedrivers, t)
    end

    b = cancercell.binitial
    d = cancercell.dinitial

    killcell = false

    #change copy number state of chromosomes
    for i in 1:cancercell.chromosomes.N
        cancercell.chromosomes.CN[i] += mutations_gain[i] - mutations_loss[i]
        #CN cannot go below minCN, cell dies
        if cancercell.chromosomes.CN[i] < minCN
            #cancercell.chromosomes.CN[i] = minCN
            killcell = true
        end
        #CN cannot exceed maxCN
        if cancercell.chromosomes.CN[i] > maxCN
            cancercell.chromosomes.CN[i] = maxCN
        end
    end

    b, d = fitnessfunc(cancercell, s, b, d)

    cancercell.b = b
    cancercell.d = d

    if cancercell.b + cancercell.d > Rmax
      Rmax = cancercell.b + cancercell.d
    end

    return cancercell, Rmax, killcell
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
    N0 = 1, μ = Chrmutrate(Nchr, m = 0.01), s = Chrfitness(Nchr, m = 0.01),
    timefunction::Function = exptime, fitnessfunc = optimumfitness(),
    maxCN = 6, minCN = 1, states = [],
    verbose = true)

    #Rmax starts with b + d and changes once a fitter mutant is introduced, this ensures that
    # b and d have correct units

    #initialize arrays and parameters
    t, tvec, N, Nvec, cells = initializesim(b, d, Nchr, N0 = N0, states = states)
    cells = getfitness(cells, s, b, d, fitnessfunc = fitnessfunc)
    brate = cells[1].b
    drate = cells[1].d
    Rmax = brate + drate
    if verbose
        println("##################################")
        println("Birth rate = $brate, death rate = $drate")
        println("initial Rmax: $Rmax")
        println("Mean fitness = $(meanfitness(cells)), Max fitness = $(maxfitness(cells)), Min fitness = $(minfitness(cells))")
        println("Initial distance from optimum: $(cells[1].chromosomes.CN .- s.optimum)")
        #println(cells[1])
        println("##################################")
    end

    while N < Nmax
        #pick a random cell
        randcell = rand(1:N)
        r = rand(Uniform(0, Rmax))
	    Nt = N
        Rmaxt = Rmax
        #println("b: $(cells[randcell].b), d: $(cells[randcell].d)")
        #println("CN: $(cells[randcell].chromosomes)")
        #println([cells[randcell].b, cells[randcell].d, r])
        #birth event if r<birthrate, access correct birthrate from cells array

        #nothing if r > b+d
        if ((cells[randcell].b + cells[randcell].d) <= r )
          push!(Nvec, N)
          Δt =  1/(Rmax * Nt) * timefunction()
          t = t + Δt
          push!(tvec,t)
        end

        #death event if b<r<b+d
        if (cells[randcell].b <= r < (cells[randcell].b + cells[randcell].d))
            #population decreases by 1
            N = N - 1
            #frequency of cell type decreases
            #remove deleted cell
            deleteat!(cells,randcell)
            push!(Nvec,N)
            Δt =  1/(Rmax * Nt) * timefunction()
            t = t + Δt
            push!(tvec,t)
            #every cell dies reinitialize simulation
            if (N == 0)
                t, tvec, N, Nvec, cells = initializesim(b, d, Nchr, N0 = N0, states = states)
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
            cells[randcell], Rmax, killcell = newmutations(cells[randcell], μ, s, Rmax,
            t, fitnessfunc, maxCN = maxCN, minCN = minCN)
            if killcell == true
                N = N - 1
                deleteat!(cells,randcell)
            end
            cells[end], Rmax, killcell = newmutations(cells[end], μ, s, Rmax, t,
            fitnessfunc, maxCN = maxCN, minCN = minCN)
            if killcell == true
                N = N - 1
                deleteat!(cells,length(cells))
            end
            push!(Nvec, N)
            Δt =  1/(Rmaxt * Nt) * timefunction()
            t = t + Δt
            push!(tvec,t)
        end

        #if randcell > length(cells)
        #    println("$N, $(length(cells)), $randcell, $r")
        #end

        #every cell dies reinitialize simulation
        if (N == 0)
            t, tvec, N, Nvec, cells = initializesim(b, d, Nchr, N0 = N0, states = states)
        end
    end
    if verbose
        println()
        println()
        println("##################################")
        println("Mean fitness = $(meanfitness(cells)), Max fitness = $(maxfitness(cells)), Min fitness = $(minfitness(cells))")
        println("Median genotype:")
        println("$(mediangenotype(cells))")
        println("Mean genotype:")
        println("$(meangenotype(cells))")
        println("Optimum genotype:")
        println("$(s.optimum)")
        println("Difference in genotype:")
        println("$(s.optimum .-mediangenotype(cells))")
        println("##################################")
        println()
    end
    return cells, (tvec, Nvec), Rmax
end



function simulate(cells::Array{cancercellCN, 1}, tvec, Nvec, Nmax;
    μ = Chrmutrate(Nchr, m = 0.01), s = Chrfitness(Nchr, m = 0.01),
    timefunction::Function = exptime, fitnessfunc = optimumfitness(),
    maxCN = 6, minCN = 1, states = [],
    verbose = true)

    #Rmax starts with b + d and changes once a fitter mutant is introduced, this ensures that
    # b and d have correct units

    t = tvec[end]
    Nvec[end] = length(cells)
    N = Nvec[end]

    #initialize arrays and parameters
    brate = maximum(map(x -> x.b, cells))
    drate = maximum(map(x -> x.d, cells))
    Rmax = brate + drate
    if verbose
        println("##################################")
        println("Current time = $t, current N = $N")
        println("Birth rate = $brate, death rate = $drate")
        println("initial Rmax: $Rmax")
        println("Mean fitness = $(meanfitness(cells)), Max fitness = $(maxfitness(cells)), Min fitness = $(minfitness(cells))")
        println("Initial distance from optimum: $(cells[1].chromosomes.CN .- s.optimum)")
        #println(cells[1])
        println("##################################")
    end

    while N < Nmax
        #pick a random cell
        randcell = rand(1:N)
        r = rand(Uniform(0, Rmax))
	    Nt = N
        Rmaxt = Rmax
        #println("b: $(cells[randcell].b), d: $(cells[randcell].d)")
        #println("CN: $(cells[randcell].chromosomes)")
        #println([cells[randcell].b, cells[randcell].d, r])
        #birth event if r<birthrate, access correct birthrate from cells array

        #nothing if r > b+d
        if ((cells[randcell].b + cells[randcell].d) <= r )
          push!(Nvec, N)
          Δt =  1/(Rmax * Nt) * timefunction()
          t = t + Δt
          push!(tvec,t)
        end

        #death event if b<r<b+d
        if (cells[randcell].b <= r < (cells[randcell].b + cells[randcell].d))
            #population decreases by 1
            N = N - 1
            #frequency of cell type decreases
            #remove deleted cell
            deleteat!(cells,randcell)
            push!(Nvec,N)
            Δt =  1/(Rmax * Nt) * timefunction()
            t = t + Δt
            push!(tvec,t)
            #every cell dies reinitialize simulation
            if (N == 0)
                t, tvec, N, Nvec, cells = initializesim(b, d, Nchr, N0 = N0, states = states)
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
            cells[randcell], Rmax, killcell = newmutations(cells[randcell], μ, s, Rmax,
            t, fitnessfunc, maxCN = maxCN, minCN = minCN)
            if killcell == true
                N = N - 1
                deleteat!(cells,randcell)
            end
            cells[end], Rmax, killcell = newmutations(cells[end], μ, s, Rmax, t,
            fitnessfunc, maxCN = maxCN, minCN = minCN)
            if killcell == true
                N = N - 1
                deleteat!(cells,length(cells))
            end
            push!(Nvec, N)
            Δt =  1/(Rmaxt * Nt) * timefunction()
            t = t + Δt
            push!(tvec,t)
        end

        #if randcell > length(cells)
        #    println("$N, $(length(cells)), $randcell, $r")
        #end

        #every cell dies reinitialize simulation
        if (N == 0)
            t, tvec, N, Nvec, cells = initializesim(b, d, Nchr, N0 = N0, states = states)
        end
    end
    if verbose
        println()
        println()
        println("##################################")
        println("Current time = $t, current N = $N")
        println("Mean fitness = $(meanfitness(cells)), Max fitness = $(maxfitness(cells)), Min fitness = $(minfitness(cells))")
        println("Median genotype:")
        println("$(mediangenotype(cells))")
        println("Mean genotype:")
        println("$(meangenotype(cells))")
        println("Optimum genotype:")
        println("$(s.optimum)")
        println("Difference in genotype:")
        println("$(s.optimum .-mediangenotype(cells))")
        println("##################################")
        println()
    end
    return cells, (tvec, Nvec), Rmax
end
