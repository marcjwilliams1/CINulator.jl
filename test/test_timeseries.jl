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

@time simresult = simulate(b, d, Nmax, Nchr,
                μ = μ, s = s,
                fitnessfunc = fopt,
                states = fill(2, Nchr))

sampledcells = samplecells(simresult.cells, 0.1)
@test length(sampledcells) == 0.1 * length(simresult.cells)

sampledcells = samplecells(simresult.cells, 100)
@test length(sampledcells) == 100

simresult = simulate(sampledcells, simresult.t, simresult.N, Nmax,
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


global simresult = simulate(b, d, Nmax, Nchr,
                μ = μ, s = s,
                fitnessfunc = fopt,
                states = fill(2, Nchr))
global sampledcells = samplecells(simresult.cells, 0.1)
for i in 1:5
    global simresult = simulate(sampledcells,
        simresult.t,
        simresult.N, Nmax,
        μ = μ, s = s,
        fitnessfunc = fopt)
    global sampledcells = samplecells(simresult.cells, 0.1)
end

Ntimepoints = 10
maxN = 10^3
cellst = simulate_timeseries(b, d, maxN, Nchr, N0 = 1,
                        Ntimepoints = Ntimepoints, s = s, μ = μ,
                        states = fill(2, Nchr), pct = 0.01)
Nvec = cellst[end].N
tvec = cellst[end].t
@test sum(abs.(Nvec[1:end-1] - Nvec[2:end]) .>100) == Ntimepoints - 1
@test sum(Nvec[1:end-1][abs.(Nvec[1:end-1] - Nvec[2:end]) .>1] .+ 1 .== maxN) == Ntimepoints - 1
plot(tvec, Nvec)


Ntimepoints = 2
simresult = simulate_timeseries(b, d, 10^5, Nchr, N0 = 1,
                        Ntimepoints = Ntimepoints, s = s, μ = μ,
                        states = fill(2, Nchr), pct = 0.1, verbose = true)
Nvec = simresult[end].N
@test sum(abs.(Nvec[1:end-1] - Nvec[2:end]) .>100) == Ntimepoints - 1
