#######################################################################################################
## Test optimum fitness model, where fitness is a function of the distance from some optimum karyotype
######################################################################################################

Random.seed!(1)

fopt = CINulator.optimumfitness()
chrfitness = CINulator.Chrfitness(2)
chrfitness.optimum .= 2 #fix high optimum so it's not reached
x = CINulator.initializesim(log(2), 0.0, 2)
cell = x[5][1]
cell.genome.CN[1] = CINulator.Chromosome(2,2)
cell.genome.CN[2] = CINulator.Chromosome(1,0)
chrfitness.optimum[1] = 4
chrfitness.optimum[2] = 1
chrfitness.alpha = fill(0.5, 2)
bdrates = fopt(cell, chrfitness, log(2), 0.1)
@test isapprox(bdrates[1], log(2), rtol = 0.01)

###################################################################################
# Test that cells with karyootypes far away from the optimum have lower birth rates
###################################################################################
cell.genome.CN[1] = CINulator.Chromosome(1,1)
cell.genome.CN[2] = CINulator.Chromosome(1,1)
bdrates1 = fopt(cell, chrfitness, log(2), 0.1)

cell.genome.CN[1] = CINulator.Chromosome(2,1)
cell.genome.CN[2] = CINulator.Chromosome(1,1)
bdrates2 = fopt(cell, chrfitness, log(2), 0.1)

cell.genome.CN[1] = CINulator.Chromosome(2,1)
cell.genome.CN[2] = CINulator.Chromosome(1,0)
bdrates3 = fopt(cell, chrfitness, log(2), 0.1)

cell.genome.CN[1] = CINulator.Chromosome(3,1)
cell.genome.CN[2] = CINulator.Chromosome(1,0)
bdrates4 = fopt(cell, chrfitness, log(2), 0.1)
@test bdrates4[1] > bdrates3[1] > bdrates2[1] > bdrates1[1]

###################################################################################
#Test full simulation using optimum
###################################################################################
Nchr = 4
Nmax = 10^5
b = log(2)
d = 0.0
μ = CINulator.Chrmutrate(Nchr, mutrates = [0.2, 0.2, 0.2, 0.2])
s = CINulator.Chrfitness(Nchr, optimum = [2, 4, 1, 3], alpha = 0.5)

t, tvec, N, Nvec, cells = CINulator.initializesim(b, d, Nchr, N0 = 1)

simresult = CINulator.simulate(b, d, Nmax, Nchr,
                μ = μ, s = s,
                fitnessfunc = fopt)

dist = map(x -> sum(abs.(gettotalcn(x) .- s.optimum)), simresult.cells)
fitness = map(x -> x.b - x.d, simresult.cells)

# fitness and distance from optimum should be negatively correlated
@test cor(dist, fitness) < 0.0
@test cor(dist, fitness) < cor(shuffle(fitness), shuffle(dist))
f = copynumberfrequency(simresult.cells)
#plot(cells)
@test findmax(f[:, 2])[2] - 1 == 4


#with alpha = 0.0 (no selection)
s = CINulator.Chrfitness(Nchr,
        optimum = [2, 4, 1, 3],
        alpha = 0.0)

simresult = CINulator.simulate(b, d, Nmax, Nchr,
                μ = μ, s = s,
                fitnessfunc = fopt)
dist = map(x -> sum(abs.(gettotalcn(x) .- s.optimum)), simresult.cells)
fitness = map(x -> x.b - x.d, simresult.cells)


@test abs(cor(dist, fitness)) < 0.001
@test sum(map(x -> x.b, simresult.cells) .== b) == Nmax
#f = copynumberfrequency(simresult.cells)
#plot(cells)


#### different alpha for different chromosomes
Random.seed!(123)
Nchr = 5
Nmax = 10^5
b = log(2)
d = 0.1
s = CINulator.Chrfitness(Nchr,
        optimum = fill(3, Nchr),
        alpha = 0.1)
s.alpha[1] = 0.0
s.alpha[2] = 5.0
s.alpha[3] = 0.5

mut = 0.1/Nchr
μ = CINulator.Chrmutrate(Nchr, mutrates = fill(mut, Nchr))

@time simresult = CINulator.simulate(b, d, Nmax, Nchr,
                μ = μ, s = s,
                fitnessfunc = fopt,
                states = fill((1,1), Nchr))
#plot(simresult.cells)
f = copynumberfrequency(simresult.cells)
@test findmax(f[:,2])[1] > findmax(f[:,3])[1]




#######################################################################################################
## Test additive fitness model
######################################################################################################

Random.seed!(1)

fadd = CINulator.additivefitness()
chrfitness = CINulator.Chrfitnessadditive(2, s = fill(2.0, 2))
chrfitness.optimum .= 2 #fix high optimum so it's not reached
x = CINulator.initializesim(log(2), 0.0, 2)
cell = x[5][1]
cell.genome.CN[1] = CINulator.Chromosome(2,2)
cell.genome.CN[2] = CINulator.Chromosome(1,0)
chrfitness.optimum[1] = 4
chrfitness.optimum[2] = 1
chrfitness.alpha = fill(0.5, 2)
bdrates = fopt(cell, chrfitness, log(2), 0.1)
@test isapprox(bdrates[1], log(2), rtol = 0.01)