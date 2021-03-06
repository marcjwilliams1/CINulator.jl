#######################################################################################################
## Test optimum fitness model, where fitness is a function of the distance from some optimum karyotype
######################################################################################################

Random.seed!(1)

fopt = CINulator.optimumfitness()
myoptgenome = Genome(2)
myoptgenome.CN[1] = Chromosome(Arm(2,2), Arm(2,2))
myoptgenome.CN[2] = Chromosome(Arm(1,0), Arm(1,0))
chrfitness = CINulator.Chrfitness(myoptgenome)
x = CINulator.initializesim(log(2), 0.0, 2)
cell = x[5][1]
cell.genome.CN[1] = CINulator.Chromosome(Arm(2,2), Arm(2,2))
cell.genome.CN[2] = CINulator.Chromosome(Arm(1,0), Arm(1,0))
chrfitness.alpha = fill(0.5, 4)
bdrates = fopt(cell, chrfitness, log(2), 0.1)
@test isapprox(bdrates[1], log(2), rtol = 0.01)

myoptgenome = Genome(22)
myoptgenome.CN[1] = Chromosome(Arm(2,2), Arm(2,2))
myoptgenome.CN[2] = Chromosome(Arm(1,0), Arm(1,0))
myoptgenome.CN[5] = Chromosome(Arm(1,1), Arm(0,2))
myoptgenome.CN[18] = Chromosome(Arm(3,1), Arm(1,0))
chrfitness = CINulator.Chrfitness(myoptgenome)
chrfitness.alpha = fill(0.1, 22*2)
CINulator.genomedistance(Genome(22), chrfitness.optimum, chrfitness.alpha)

###################################################################################
# Test that cells with karyootypes far away from the optimum have lower birth rates
###################################################################################
fopt = CINulator.optimumfitness()
myoptgenome = Genome(2)
myoptgenome.CN[1] = Chromosome(Arm(2,2), Arm(2,2))
myoptgenome.CN[2] = Chromosome(Arm(1,0), Arm(1,0))
chrfitness = CINulator.Chrfitness(myoptgenome)
chrfitness.alpha = fill(0.5, 4)

cell.genome.CN[1] = CINulator.Chromosome(Arm(1,1), Arm(1,1))
cell.genome.CN[2] = CINulator.Chromosome(Arm(1,1), Arm(1,1))
bdrates1 = fopt(cell, chrfitness, log(2), 0.1)

cell.genome.CN[1] = CINulator.Chromosome(Arm(2,1), Arm(2,1))
cell.genome.CN[2] = CINulator.Chromosome(Arm(1,1), Arm(1,1))
bdrates2 = fopt(cell, chrfitness, log(2), 0.1)

cell.genome.CN[1] = CINulator.Chromosome(Arm(2,1), Arm(2,1))
cell.genome.CN[2] = CINulator.Chromosome(Arm(1,0), Arm(1,0))
bdrates3 = fopt(cell, chrfitness, log(2), 0.1)

cell.genome.CN[1] = CINulator.Chromosome(Arm(3,1), Arm(3,1))
cell.genome.CN[2] = CINulator.Chromosome(Arm(1,0), Arm(1,0))
bdrates4 = fopt(cell, chrfitness, log(2), 0.1)
@test bdrates4[1] > bdrates3[1] > bdrates2[1] > bdrates1[1]

###################################################################################
# Test when chromosome have different alpha
###################################################################################
fopt = CINulator.optimumfitness()
foptmod = CINulator.optimumfitness(distancemethod = "mod")
myoptgenome = Genome(2)
myoptgenome.CN[1] = Chromosome(Arm(2,2), Arm(2,2))
myoptgenome.CN[2] = Chromosome(Arm(1,0), Arm(1,0))
chrfitness = CINulator.Chrfitness(myoptgenome, alpha = [0.5, 0.0, 0.5, 0.0])

#Same karyotype, no change in fitness
cell.genome.CN[1] = CINulator.Chromosome(Arm(2,2), Arm(2,2))
cell.genome.CN[2] = CINulator.Chromosome(Arm(1,0), Arm(1,0))
bdrates1 = fopt(cell, chrfitness, log(2), 0.1)
bdrates1mod = foptmod(cell, chrfitness, log(2), 0.1)

#Different karyotype + no change in fitness
cell.genome.CN[1] = CINulator.Chromosome(Arm(2,2), Arm(1,1))
cell.genome.CN[2] = CINulator.Chromosome(Arm(1,0), Arm(1,1))
bdrates2 = fopt(cell, chrfitness, log(2), 0.1)
bdrates2mod = foptmod(cell, chrfitness, log(2), 0.1)

#Different karyotype + change in fitness
cell.genome.CN[1] = CINulator.Chromosome(Arm(2,1), Arm(1,1))
cell.genome.CN[2] = CINulator.Chromosome(Arm(1,0), Arm(1,1))
bdrates3 = fopt(cell, chrfitness, log(2), 0.1)
bdrates3mod = foptmod(cell, chrfitness, log(2), 0.1)

@test bdrates3[1] != bdrates3mod[1]
@test bdrates1[1] == bdrates2[1]
@test bdrates3[1] < bdrates2[1]

###################################################################################
#Test full simulation using optimum
###################################################################################
Nchr = 4
Nmax = 10^5
b = log(2)
d = 0.0

myoptgenome = Genome(4)
myoptgenome.CN[2] = Chromosome(Arm(2,2), Arm(2,2))
myoptgenome.CN[3] = Chromosome(Arm(1,0), Arm(1,0))
myoptgenome.CN[4] = Chromosome(Arm(2,1), Arm(2,1))
s = CINulator.Chrfitness(myoptgenome, alpha = 0.75)
μ = CINulator.Chrmutrate(Nchr, mutrates = fill(0.2, Nchr), psplit = fill(0.0, Nchr))
fopt = CINulator.optimumfitness()

t, tvec, N, Nvec, cells = CINulator.initializesim(b, d, Nchr, N0 = 1)

simresult = CINulator.simulate(b, d, Nmax, Nchr,
                μ = μ, s = s,
                fitnessfunc = fopt, states = Genome(Nchr))

dist = map(x -> sum(abs.(gettotalcn(x) .- gettotalcn(s.optimum))), simresult.cells)
fitness = map(x -> x.b - x.d, simresult.cells)

# fitness and distance from optimum should be negatively correlated
@test cor(dist, fitness) < 0.0
@test cor(dist, fitness) < cor(shuffle(fitness), shuffle(dist))
f = copynumberfrequency(simresult.cells)
#plot(cells)
@test findmax(f[:, 3])[2] - 1 == 4
@test f[:,1] == f[:,2]

# check all cells have positive copy number for all chromosomes
df = mergecelldataframe_locus(simresult.cells)
@test sum(df[!, :state] .< 1) == 0


#with alpha = 0.0 (no selection)
s = CINulator.Chrfitness(myoptgenome, alpha = 0.0)

simresult = CINulator.simulate(b, d, Nmax, Nchr,
                μ = μ, s = s,
                fitnessfunc = fopt)
dist = map(x -> sum(abs.(gettotalcn(x) .- gettotalcn(s.optimum))), simresult.cells)
fitness = map(x -> x.b - x.d, simresult.cells)


@test abs(cor(dist, fitness)) < 0.001
@test sum(map(x -> x.b, simresult.cells) .== b) == Nmax
#f = copynumberfrequency(simresult.cells)
#plot(cells)

###################################################################################
#Include split chromosomes
###################################################################################
μ = CINulator.Chrmutrate(Nchr, mutrates = fill(0.2, Nchr), psplit = fill(0.1, Nchr))
s = CINulator.Chrfitness(myoptgenome, alpha = 0.5)
simresult = CINulator.simulate(b, d, Nmax, Nchr,
                μ = μ, s = s,
                fitnessfunc = fopt)
f = copynumberfrequency(simresult.cells)
@test f[:,1] != f[:,2]
@test findmax(f[:, 3])[2] - 1 == 4
@test findmax(f[:, 4])[2] - 1 == 4


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

#######################################################################################################
## 22 chromosomes
######################################################################################################

Random.seed!(3)
Nchr = 22
Nmax = 10^3
b = log(2)
d = 0.0
mu = 0.1

myoptgenome = Genome(Nchr)
s = CINulator.Chrfitness(myoptgenome, alpha = 0.0001)
μ = CINulator.Chrmutrate(Nchr, mutrates = fill(mu, Nchr), psplit = fill(0.1, Nchr))
fopt = CINulator.optimumfitness()
Random.seed!(3)
simresult = CINulator.simulate(b, d, Nmax, Nchr,
                μ = μ, s = s,
                fitnessfunc = fopt, states = myoptgenome)
df = mergecelldataframe_locus(simresult.cells)
@linq by(df,:cell_id, Aloh = sum(:A .== 0), Bloh = sum(:B .== 0)) |> orderby(:Bloh)
df2 = deepcopy(df)
using DataFrames
using DataFramesMeta
dfsummary = @linq by(df,:cell_id, Aloh = sum(:A .== 0), Bloh = sum(:B .== 0)) |> orderby(:Bloh)

@where(df, :cell_id .== dfsummary[!,:cell_id][end])
#|> where(:Aloh .> 0 | :Bloh .> 0)





