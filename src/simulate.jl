mutable struct Chromosomes
    CN::Array{Int64, 1}
    N::Int64
    function Chromosomes(N)
        return new(fill(2, N), N)
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
        if isempty(mutratesgain)
            mutratesloss = rand(Gamma(m), Int64(N))
            mutratesgain = rand(Gamma(m), Int64(N))
        end
        return new(mutratesloss, mutratesgain)
    end
end

mutable struct Chrfitness
    loss::Array{Float64, 1}
    gain::Array{Float64, 1}
    function Chrfitness(N; m = 0.1, fitnessloss = [], fitnessgain = [])
        if isempty(fitnessloss)
            fitnessloss = rand(Gamma(m), N)
            fitnessgain = rand(Gamma(m), N)
        end
        return new(fitnessloss, fitnessgain)
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

function initializesim(b, d, Nchr)

    #initialize time to zero
    t = 0.0
    tvec = Float64[]
    push!(tvec,t)

    #population starts with one cell
    N = 1
    Nvec = Int64[]
    push!(Nvec,N)

    #Initialize array of cell type that stores mutations for each cell and their fitness type
    #fitness type of 1 is the host population, lowest fitness
    cells = cancercellCN[]
    push!(cells, cancercellCN(b, d, b, d, Float64[], [], Chromosomes(Nchr)))

    return t, tvec, N, Nvec, cells
end

function multiplicativefitness(cancercell, chrfitness, b)

  gain_or_loss = cancercell.chromosomes.CN -
        fill(2, cancercell.chromosomes.N)

  j = 1
  for CNstate in gain_or_loss
      gain = CNstate > 0
      for i in 1:abs(CNstate)
          if gain == true
              b = b * (1.0 + chrfitness.gain[j])
          else
              b = b * (1.0 + chrfitness.loss[j])
          end
      end
      j += 1
  end

  return b
end


function newmutations(cancercell::cancercellCN,
    μ::Chrmutrate,
    s::Chrfitness,
    Rmax, t,
    fitnessfunc)

    #function to add new mutations to cells based on μ
    mutations_gain = map(x -> rand(Poisson(x)), μ.gain) .> 0
    mutations_loss = map(x -> rand(Poisson(x)), μ.loss) .> 0

    #record time if mutations occur
    if (sum(mutations_gain) + sum(mutations_loss)) > 0
        push!(cancercell.timedrivers, t)
    end

    b = cancercell.binitial

    #change copy number state of chromosomes
    for i in 1:cancercell.chromosomes.N
        cancercell.chromosomes.CN[i] += mutations_gain[i] - mutations_loss[i]
        if cancercell.chromosomes.CN[i] < 0
            cancercell.chromosomes.CN[i] = 0
        end
    end

    b = fitnessfunc(cancercell, s, b)

    cancercell.b = b

    if cancercell.b + cancercell.d > Rmax
      Rmax = cancercell.b + cancercell.d
    end

    return cancercell, Rmax
end

exptime() = - log(rand())
meantime() = 1

function simulate(b, d, Nmax, Nchr;
    μ = Chrmutrate(Nchr, m = 0.01), s = Chrfitness(Nchr, m = 0.01),
    timefunction::Function = exptime, fitnessfunc = multiplicativefitness)

    #Rmax starts with b + d and changes once a fitter mutant is introduced, this ensures that
    # b and d have correct units
    Rmax = b + d

    #initialize arrays and parameters
    t, tvec, N, Nvec, cells = initializesim(b, d, Nchr)
    while N < Nmax

        #pick a random cell
        randcell = rand(1:N)
        r = rand(Uniform(0, Rmax))
	    Nt = N
        Rmaxt = Rmax

        #birth event if r<birthrate, access correct birthrate from cells array
        if r < cells[randcell].b

            #population increases by one
            N = N + 1
            #copy cell and mutations for cell that reproduces
            push!(cells,copycell(cells[randcell]))
            #add new mutations to both new cells
            cells[randcell], Rmax = newmutations(cells[randcell], μ, s, Rmax,
            t, fitnessfunc)
            cells[end], Rmax = newmutations(cells[end], μ, s, Rmax, t,
            fitnessfunc)
            push!(Nvec, N)
            Δt =  1/(Rmaxt * Nt) * timefunction()
            t = t + Δt
            push!(tvec,t)
        end

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
        end

        #every cell dies reinitialize simulation
        if (N == 0)
            t, tvec, N, Nvec, cells = initializesim(b, d, Nchr)
        end
    end

    return cells, tvec, Rmax
end

function copynumberfrequency(cells::Array{cancercellCN, 1})

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