Random.seed!(1)
fopt = CINulator.optimumfitness()
Nchr = 5
Nmax = 10^5
b = log(2)
d = 0.0
s = CINulator.Chrfitness(Nchr,
        optimum = fill(3, Nchr),
        alpha = 0.0)
s.alpha[1] = 0.0
s.alpha[2] = 0.0
s.alpha[3] = 0.0

mut = 0.1/Nchr
μ = CINulator.Chrmutrate(Nchr,
        mutratesgain = fill(mut, Nchr),
        mutratesloss = fill(mut, Nchr))

@time simresult = simulate(b, d, 2^6, Nchr,
                μ = μ, s = s,
                fitnessfunc = fopt,
                states = fill(2, Nchr),
                timestop = true,
                timefunction = CINulator.exptime)

global Nvec = Int64[]
global tvec = Float64[]

for i in 1:1000
    @time simresult = simulate(b, d, 2^6, Nchr,
                    μ = μ, s = s,
                    fitnessfunc = fopt,
                    states = fill(2, Nchr),
                    timestop = true, verbose = false,
                    timefunction = CINulator.exptime)
    push!(Nvec, simresult.N[end])
    push!(tvec, simresult.t[end])
end



Ntimepoints = 10
simresult = simulate_timeseries(b, d, 10^3,
                        Nchr, N0 = 1,
                        Ntimepoints = Ntimepoints,
                        s = s, μ = μ,
                        tend = 0.0,
                        timestop = true,
                        timefunction = CINulator.exptime,
                        states = fill(2, Nchr),
                        pct = 0.1)
@test sum(abs.(simresult[end].N[1:end-1] - simresult[end].N[2:end]) .>1) == Ntimepoints





################################
# Approx fixed population
################################

Random.seed!(1)
fopt = CINulator.optimumfitness()
Nchr = 5
N0 = 10^3
Nmax = 10^5
b = log(2)
d = log(2)
s = CINulator.Chrfitness(Nchr,
        optimum = fill(3, Nchr),
        alpha = 0.0)
s.alpha[1] = 0.0
s.alpha[2] = 0.0
s.alpha[3] = 0.0

mut = 0.4/Nchr
μ = CINulator.Chrmutrate(Nchr,
        mutratesgain = fill(mut, Nchr),
        mutratesloss = fill(mut, Nchr))

@time simresult = simulate(b, d, Nmax, Nchr,
                μ = μ, s = s, N0 = N0,
                fitnessfunc = fopt,
                states = fill(2, Nchr),
                timestop = true,
                tend = 5.0,
                maxCN = 8,
                record = true,
                timefunction = CINulator.exptime)
plot(simresult.t, simresult.fitness)
