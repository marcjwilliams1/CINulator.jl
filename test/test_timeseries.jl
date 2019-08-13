Random.seed!(1)
fopt = CINulator.optimumfitness()
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

sampledcells = samplecells(cells, 0.1)
@time length(sampledcells) == 0.1 * length(cells)

sampledcells = samplecells(cells, 100)
@time length(sampledcells) == 100
