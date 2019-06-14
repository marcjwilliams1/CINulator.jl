Random.seed!(1)

x = CINulator.Chromosomes(20)
@test length(x.CN) == 20

mu = CINulator.Chrmutrate(20, m = 0.1)
@test length(mu.gain) == 20
@test (mu.gain .>= 0.0) == fill(true, length(mu.gain))

fit = CINulator.Chrfitness(100, m = 0.1)
@test length(fit.gain) == 100
@test (fit.gain .>= 0.0) == fill(true, length(fit.gain))
@test (fit.loss .< fit.gain) == (fit.optimum .!= 1)

x = CINulator.initializesim(0.5, 0.0, 20)
fit = CINulator.Chrfitness(20, m = 0.1)
fit.optimum .= 10 #fix high optimum so it's not reached
@test x[5][1].chromosomes.CN[1] == 2

#test fitness calculations
x[5][1].chromosomes.CN[1] = 3
x[5][1].chromosomes.CN[3] = 1
fit.gain[1] = 0.5
fit.loss[3] = -0.9
bdrates = CINulator.multiplicativefitness(x[5][1], fit, log(2), 0.1)
@test isapprox(bdrates[1], (log(2) * (1 + 0.5)), rtol = 0.01)
@test isapprox(bdrates[2], 0.1 * (1 + 0.9))

x = CINulator.initializesim(0.5, 0.0, 20)
x[5][1].chromosomes.CN[1] = 4
x[5][1].chromosomes.CN[3] = 1
fit.gain[1] = 0.5
fit.loss[3] = -0.9
bdrates = CINulator.multiplicativefitness(x[5][1], fit, log(2), 0.1)
@test isapprox(bdrates[1], (log(2) * (1 + 0.5)^2), rtol = 0.01)
@test isapprox(bdrates[2], (0.1 * (1 +  0.9)), rtol = 0.01)

x = CINulator.initializesim(0.5, 0.0, 20)
x[5][1].chromosomes.CN[1] = 3
x[5][1].chromosomes.CN[3] = 0
fit.gain[1] = 0.5
fit.loss[3] = -0.9
bdrates = CINulator.multiplicativefitness(x[5][1], fit, log(2), 0.1)
@test isapprox(bdrates[1], (log(2) * (1 + 0.5)^1), rtol = 0.01)
@test isapprox(bdrates[2], (0.1 * (1 +  0.9) ^2), rtol = 0.01)

x = CINulator.initializesim(0.5, 0.0, 20)
x[5][1].chromosomes.CN[1] = 2
x[5][1].chromosomes.CN[3] = 2
mu.gain[3] = 0.0
mu.loss[1] = 0.0
mu.gain[1] = 10.0
mu.loss[3] = 10.0
c, Rmax = CINulator.newmutations(x[5][1], mu, fit, 1.0, 10.0, CINulator.multiplicativefitness)
@test c.chromosomes.CN[1] == 3
@test c.chromosomes.CN[3] == 1

#test if passing optimum results in decreased fitness
x = CINulator.initializesim(0.5, 0.0, 20)
fit = CINulator.Chrfitness(20, m = 0.1)
fit.optimum .= 4

x[5][1].chromosomes.CN[1] = 5
x[5][1].chromosomes.CN[3] = 1
fit.gain[1] = 0.5
fit.loss[3] = -0.9
bdrates = CINulator.multiplicativefitness(x[5][1], fit, log(2), 0.1)
@test isapprox(bdrates[1], (log(2) * (1 + 0.5)^3 * (1 + 0.5) ^ -2), rtol = 0.01)
@test isapprox(bdrates[2], 0.1 * (1 + 0.9))

#copy number state of 3 should have same fitness as copy number state 5 because
# the optimum is 4
x[5][1].chromosomes.CN[1] = 3
bdrates2 = CINulator.multiplicativefitness(x[5][1], fit, log(2), 0.1)
@test isapprox(bdrates2[1], (log(2) * (1 + 0.5)^1 * (1 + 0.5) ^ 0), rtol = 0.01)
@test bdrates2[1] == bdrates[1]

#
rates = Float64[]
for i in [2,3,4,5,6,7]
    x[5][1].chromosomes.CN[1] = i
    bdrates = CINulator.multiplicativefitness(x[5][1], fit, log(2), 0.1)
    push!(rates, bdrates[1])
end




#test simulation
Nchr = 4
Nmax = 10^5
b = log(2)
d = 0.2
μ = CINulator.Chrmutrate(Nchr,
        mutratesgain = [0.2, 0.2, 0.2, 0.0],
        mutratesloss = [0.2, 0.2, 0.2, 0.0])
s = CINulator.Chrfitness(Nchr, fitnessgain = [0.0, 0.1, 0.0, 0.0],
        fitnessloss = [0.0, 0.0, 0.1, 0.0],
        optimum = [10, 10, 10, 10])
cells, t, Rmax = CINulator.simulate(b, d, Nmax, Nchr,
                μ = μ, s = s)
f = copynumberfrequency(cells)
plot(cells)

#Chromosome 4 has 0 mutation rate so all cells should have CN = 2
@test f[3, Nchr] == 1.0

#average copy  number in chr3 < chr1 < chr2
@test sum(f[:, 3] .* collect(0:size(f)[1] - 1)) <
    sum(f[:, 1] .* collect(0:size(f)[1] - 1)) <
    sum(f[:, 2] .* collect(0:size(f)[1] - 1))

#Chromosome 2, gains are positively selected, modal value of CN should be > 2
@test findmax(f[:, 2])[2] - 1 > 2

#Chromosome 3, losses are positively selected, modal value of CN should be < 2
@test findmax(f[:, 3])[2] - 1 < 2



#test simulation
Nchr = 4
Nmax = 10^5
b = log(2)
d = 0.3
μ = CINulator.Chrmutrate(Nchr,
        mutratesgain = [0.2, 0.2, 0.2, 0.0],
        mutratesloss = [0.2, 0.2, 0.2, 0.0])
s = CINulator.Chrfitness(Nchr, fitnessgain = [0.0, -0.9, 0.0, 0.0],
        fitnessloss = [0.0, -0.9, 0.1, 0.0],
        optimum = [10, 10, 10, 10])
cells, t, Rmax = CINulator.simulate(b, d, Nmax, Nchr,
                μ = μ, s = s)
f = copynumberfrequency(cells)
plot(cells)

#imposed strong selection for chromosome 2 to be diploid, so modal value should be 2.
@test findmax(f[:, 2])[2] - 1 == 2

#test adding optimum
Nchr = 4
Nmax = 10^5
b = log(2)
d = 0.2
μ = CINulator.Chrmutrate(Nchr,
        mutratesgain = [0.2, 0.2, 0.2, 0.0],
        mutratesloss = [0.2, 0.2, 0.2, 0.0])
s = CINulator.Chrfitness(Nchr, fitnessgain = [0.1, 0.3, 0.0, 0.0],
        fitnessloss = [0.1, 0.0, 0.3, 0.0],
        optimum = [2, 4, 1, 10])
cells, t, Rmax = CINulator.simulate(b, d, Nmax, Nchr,
                μ = μ, s = s)
f = copynumberfrequency(cells)
plot(cells)
@test findmax(f[:, 2])[2] - 1 == 4



# test if initializing with N0 > 1 works
N0 = 100
t, tvec, N, Nvec, cells = CINulator.initializesim(log(2), 0.0, 5, N0 = N0)
@test length(cells) == N0
cells1, t1, Rmax1 = CINulator.simulate(b, d, 200, Nchr,
                μ = μ, s = s, N0 = 1)
cells2, t2, Rmax2 = CINulator.simulate(b, d, 200, Nchr,
                μ = μ, s = s, N0 = 100)
@test length(t2) < length(t1)
