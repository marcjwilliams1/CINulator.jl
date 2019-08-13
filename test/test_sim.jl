Random.seed!(1)

x = CINulator.Chromosomes(20)
@test length(x.CN) == 20

mu = CINulator.Chrmutrate(20)
@test length(mu.gain) == 20
@test (mu.gain .>= 0.0) == fill(true, length(mu.gain))

s = CINulator.Chrfitness(100)
@test length(s.optimum) == 100

x = CINulator.initializesim(0.5, 0.0, 2)
s = CINulator.Chrfitness(2)
s.optimum .= 2 #fix high optimum so it's not reached
@test x[5][1].chromosomes.CN[1] == 2

fopt = CINulator.optimumfitness()
x = CINulator.initializesim(0.5, 0.0, 2)
x[5][1].chromosomes.CN[1] = 2
x[5][1].chromosomes.CN[2] = 2
mu.gain[2] = 0.0
mu.loss[1] = 0.0
mu.gain[1] = 10.0
mu.loss[2] = 10.0
c, Rmax = CINulator.newmutations(x[5][1], mu, s, 1.0, 10.0, fopt)
@test c.chromosomes.CN[1] == 3
@test c.chromosomes.CN[2] == 1

#test if passing optimum results in decreased fitness
x = CINulator.initializesim(0.5, 0.0, 2)
s = CINulator.Chrfitness(2)
s.optimum .= 4
