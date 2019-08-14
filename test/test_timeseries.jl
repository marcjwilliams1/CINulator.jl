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

@time cellst1, t, Rmax = simulate(b, d, Nmax, Nchr,
                μ = μ, s = s,
                fitnessfunc = fopt,
                states = fill(2, Nchr))

sampledcells = samplecells(cellst1, 0.1)
@time length(sampledcells) == 0.1 * length(cells)

sampledcells = samplecells(cells, 100)
@time length(sampledcells) == 100

cellst2, t, Rmax = simulate(sampledcells, t[1], t[2], Nmax,
                        μ = μ, s = s,
                        fitnessfunc = fopt)


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
s.alpha[2] = 0.5
s.alpha[3] = 0.5
Nmax = 10^2
mut = 0.1/Nchr
μ = CINulator.Chrmutrate(Nchr,
        mutratesgain = fill(mut, Nchr),
        mutratesloss = fill(mut, Nchr))


global cellst1, t, Rmax = simulate(b, d, Nmax, Nchr,
                μ = μ, s = s,
                fitnessfunc = fopt,
                states = fill(2, Nchr))
global sampledcells = samplecells(cellst1, 0.1)
for i in 1:5
    global cellst2, t2, Rmax = simulate(sampledcells, t[1],
        t[2], Nmax,
        μ = μ, s = s,
        fitnessfunc = fopt)
    global sampledcells = samplecells(cellst2, 0.1)
end

Ntimepoints = 10
maxN = 10^3
cellst, (tvec, Nvec) = simulate_timeseries(b, d, maxN, Nchr, N0 = 1,
                        Ntimepoints = Ntimepoints, s = s, μ = μ,
                        states = fill(2, Nchr), pct = 0.01)
@test sum(abs.(Nvec[1:end-1] - Nvec[2:end]) .>1) == Ntimepoints
@test sum(Nvec[1:end-1][abs.(Nvec[1:end-1] - Nvec[2:end]) .>1] .+ 1 .== maxN) == Ntimepoints


Ntimepoints = 2
cellst, (tvec, Nvec) = simulate_timeseries(b, d, 10^6, Nchr, N0 = 1,
                        Ntimepoints = Ntimepoints, s = s, μ = μ,
                        states = fill(2, Nchr), pct = 0.1, verbose = false)
@test sum(abs.(Nvec[1:end-1] - Nvec[2:end]) .>1) == Ntimepoints
