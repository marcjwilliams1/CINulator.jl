Nchr = 10
mut = 0.1/Nchr
μ = CINulator.Chrmutrate(Nchr,
        mutratesgain = fill(mut, Nchr),
        mutratesloss = fill(mut, Nchr))
b = log(2)
d = 0.1
s = CINulator.Chrfitness(Nchr,
        optimum = fill(3, Nchr),
        alpha = 0.1)
s.alpha[1] = 0.0
s.alpha[2] = 5.0
s.alpha[3] = 0.5
Nmax = 10^3
fopt = CINulator.optimumfitness()

@time simresult = CINulator.simulate(b, d, Nmax, Nchr,
                μ = μ, s = s,
                fitnessfunc = fopt,
                states = fill(2, Nchr),
                labelcells = true)
