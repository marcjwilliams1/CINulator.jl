Random.seed!(1)

fopt = CINulator.optimumfitness()

chrfitness = CINulator.Chrfitness(2)
chrfitness.optimum .= 2 #fix high optimum so it's not reached
x = CINulator.initializesim(log(2), 0.0, 2)
x[5][1].chromosomes.CN[1] = 4
x[5][1].chromosomes.CN[2] = 1
chrfitness.optimum[1] = 4
chrfitness.optimum[2] = 1
chrfitness.alpha = fill(0.5, 2)
bdrates = fopt(x[5][1], chrfitness, log(2), 0.1)
@test isapprox(bdrates[1], log(2), rtol = 0.01)


x[5][1].chromosomes.CN[1] = 2
x[5][1].chromosomes.CN[2] = 2
bdrates1 = fopt(x[5][1], chrfitness, log(2), 0.1)

x[5][1].chromosomes.CN[1] = 3
x[5][1].chromosomes.CN[2] = 2
bdrates2 = fopt(x[5][1], chrfitness, log(2), 0.1)

x[5][1].chromosomes.CN[1] = 3
x[5][1].chromosomes.CN[2] = 1
bdrates3 = fopt(x[5][1], chrfitness, log(2), 0.1)

x[5][1].chromosomes.CN[1] = 4
x[5][1].chromosomes.CN[2] = 1
bdrates4 = fopt(x[5][1], chrfitness, log(2), 0.1)

@test bdrates4[1] > bdrates3[1] > bdrates2[1] > bdrates1[1]


#test adding optimum
Nchr = 4
Nmax = 10^5
b = log(2)
d = 0.0
μ = CINulator.Chrmutrate(Nchr,
        mutratesgain = [0.2, 0.2, 0.2, 0.2],
        mutratesloss = [0.2, 0.2, 0.2, 0.2])
s = CINulator.Chrfitness(Nchr,
        optimum = [2, 4, 1, 3],
        alpha = 0.5)


t, tvec, N, Nvec, cells = CINulator.initializesim(b, d, Nchr, N0 = 1)

cells, t, Rmax = CINulator.simulate(b, d, Nmax, Nchr,
                μ = μ, s = s,
                fitnessfunc = fopt)

dist = map(x -> sum(abs.(x.chromosomes.CN .- s.optimum)), cells)
fitness = map(x -> x.b - x.d, cells)

# fitness and distance from optimum should be negatively correlated
@test cor(dist, fitness) < 0.0
@test cor(dist, fitness) < cor(shuffle(fitness), shuffle(dist))
f = copynumberfrequency(cells)
plot(cells)
@test findmax(f[:, 2])[2] - 1 == 4




#with alpha = 0.0
s = CINulator.Chrfitness(Nchr,
        optimum = [2, 4, 1, 3],
        alpha = 0.0)

cells, t, Rmax = CINulator.simulate(b, d, Nmax, Nchr,
                μ = μ, s = s,
                fitnessfunc = fopt)
dist = map(x -> sum(abs.(x.chromosomes.CN .- s.optimum)), cells)
fitness = map(x -> x.b - x.d, cells)


@test abs(cor(dist, fitness)) < 0.001
@test sum(map(x -> x.b, cells) .== b) == Nmax
f = copynumberfrequency(cells)
plot(cells)




#with alpha = 0.0, small Nmax
Nmax = 10^2
Nchr = 4
s = CINulator.Chrfitness(Nchr,
        optimum = [2, 2, 2, 2],
        alpha = 1.0)
mut = 0.01
μ = CINulator.Chrmutrate(Nchr,
        mutratesgain = fill(mut, 4),
        mutratesloss = fill(mut, 4))
d = 0.1

cells, t, Rmax = CINulator.simulate(b, d, Nmax, Nchr,
                μ = μ, s = s,
                fitnessfunc = fopt)
@test abs(cor(dist, fitness)) < 0.001
f = copynumberfrequency(cells)
plot(cells)

df = celldataframe(cells[1])
df = mergecelldataframe(cells)

df1 = CINulator.copynumbernoise(df, celldist = Beta(0.01, 0.1))
df2 = CINulator.copynumbernoise(df, celldist = Beta(1, 0.1))

#check initialization
mut = 0.00
μ = CINulator.Chrmutrate(Nchr,
        mutratesgain = fill(mut, 4),
        mutratesloss = fill(mut, 4))
cells, t, Rmax = CINulator.simulate(b, d, Nmax, Nchr,
                μ = μ, s = s,
                fitnessfunc = fopt,
                states = [4,3,2,1])
f = copynumberfrequency(cells)
@test findmax(f[:, 1])[2] - 1 == 4
@test findmax(f[:, 2])[2] - 1 == 3
@test findmax(f[:, 3])[2] - 1 == 2
@test findmax(f[:, 4])[2] - 1 == 1


opt = 1.0
Nchr = 22
Nmax = 10^4
b = log(2)
d = 0.3
s = CINulator.Chrfitness(Nchr,
        optimum = fill(2, Nchr),
        alpha = opt)

mut = 0.1/Nchr
μ = CINulator.Chrmutrate(Nchr,
        mutratesgain = fill(mut, Nchr),
        mutratesloss = fill(mut, Nchr))

Random.seed!(1234)
states = fill(2, Nchr)
states[1] = 1
t, tvec, N, Nvec, cells = CINulator.initializesim(b, d, Nchr, N0 = 1, states = states)
cells = CINulator.getfitness(cells, s, b, d, fitnessfunc = fopt)
@time cells, t, Rmax = CINulator.simulate(b, d, Nmax, Nchr,
                μ = μ, s = s,
                fitnessfunc = fopt)



Nchr = 22
Nmax = 10^4
b = log(2)
d = 0.1
s = CINulator.Chrfitness(Nchr,
        optimum = fill(2, Nchr),
        alpha = opt)

mut = 0.1/Nchr
μ = CINulator.Chrmutrate(Nchr,
        mutratesgain = fill(mut, Nchr),
        mutratesloss = fill(mut, Nchr))

@time cells, t, Rmax = CINulator.simulate(b, d, Nmax, Nchr,
                μ = μ, s = s,
                fitnessfunc = fopt)





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
μ = CINulator.Chrmutrate(Nchr,
        mutratesgain = fill(mut, Nchr),
        mutratesloss = fill(mut, Nchr))

@time cells, t, Rmax = CINulator.simulate(b, d, Nmax, Nchr,
                μ = μ, s = s,
                fitnessfunc = fopt,
                states = fill(2, Nchr))
plot(cells)
f = copynumberfrequency(cells)
@test findmax(f[:,2])[1] > findmax(f[:,3])[1]
