Random.seed!(1)

x = CINulator.Chromosomes(20)
@test length(x.CN) == 20

mu = CINulator.Chrmutrate(20, m = 0.1)
@test length(mu.gain) == 20
@test (mu.gain .>= 0.0) == fill(true, length(mu.gain))

fit = CINulator.Chrfitness(20, m = 0.1)
@test length(fit.gain) == 20
@test (fit.gain .>= 0.0) == fill(true, length(fit.gain))

x = CINulator.initializesim(0.5, 0.0, 20)
@test x[5][1].chromosomes.CN[1] == 2

#test fitness calculations
x[5][1].chromosomes.CN[1] = 3
x[5][1].chromosomes.CN[3] = 1
fit.gain[1] = 0.5
fit.loss[3] = -0.9
b = CINulator.multiplicativefitness(x[5][1], fit, log(2))
@test isapprox(b, (log(2) * (1 + 0.5) * (1 -  0.9)), rtol = 0.01)

x[5][1].chromosomes.CN[1] = 4
x[5][1].chromosomes.CN[3] = 1
fit.gain[1] = 0.5
fit.loss[3] = -0.9
b = CINulator.multiplicativefitness(x[5][1], fit, log(2))
@test isapprox(b, (log(2) * (1 + 0.5)^2 * (1 -  0.9)), rtol = 0.01)

x[5][1].chromosomes.CN[1] = 3
x[5][1].chromosomes.CN[3] = 0
fit.gain[1] = 0.5
fit.loss[3] = -0.9
b = CINulator.multiplicativefitness(x[5][1], fit, log(2))
@test isapprox(b, (log(2) * (1 + 0.5)^1 * (1 -  0.9) ^2), rtol = 0.01)


x[5][1].chromosomes.CN[1] = 2
x[5][1].chromosomes.CN[3] = 2
mu.gain[3] = 0.0
mu.loss[1] = 0.0
mu.gain[1] = 10.0
mu.loss[3] = 10.0
c, Rmax = CINulator.newmutations(x[5][1], mu, fit, 1.0, 10.0, CINulator.multiplicativefitness)
@test c.chromosomes.CN[1] == 3
@test c.chromosomes.CN[3] == 1


#test simulation
Nchr = 4
Nmax = 10^5
b = log(2)
d = 0.0
μ = CINulator.Chrmutrate(Nchr,
        mutratesgain = [0.2, 0.2, 0.2, 0.0],
        mutratesloss = [0.2, 0.2, 0.2, 0.0])
s = CINulator.Chrfitness(Nchr, fitnessgain = [0.0, 0.1, 0.0, 0.0],
        fitnessloss = [0.0, 0.0, 0.1, 0.0])
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
d = 0.0
μ = CINulator.Chrmutrate(Nchr,
        mutratesgain = [0.2, 0.2, 0.2, 0.0],
        mutratesloss = [0.2, 0.2, 0.2, 0.0])
s = CINulator.Chrfitness(Nchr, fitnessgain = [0.0, -0.9, 0.0, 0.0],
        fitnessloss = [0.0, -0.9, 0.1, 0.0])
cells, t, Rmax = CINulator.simulate(b, d, Nmax, Nchr,
                μ = μ, s = s)
f = copynumberfrequency(cells)
plot(cells)

#imposed strong selection for chromosome 2 to be diploid, so modal value should be 2.
@test findmax(f[:, 2])[2] - 1 == 2
