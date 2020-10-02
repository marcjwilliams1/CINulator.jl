Nchr = 2
mut = 1
μ = CINulator.Chrmutrate(Nchr, mutrates = fill(mut, Nchr))
s = CINulator.Chrfitness(Nchr,
        optimum = [2, 4],
        alpha = 0.5)
b = log(2)
d = 0.1

cell = CINulator.cancercellCN(log(2), 0.0, log(2), 0.0, Float64[], [], CINulator.Genome(Nchr, []), Int64[], "temp")
cancercell1 = CINulator.copycell(cell)
cancercell2 = CINulator.copycell(cell)

cancercell1A, cancercell2A = CINulator.missegregation(cancercell1, cancercell2, 1)
@test (map(x -> x.tot, cancercell1A.genome.CN) .+ map(x -> x.tot, cancercell2A.genome.CN)) / 2 == fill(2.0, Nchr)
@test cancercell1A.genome.CN[1] != cancercell2A.genome.CN[1]
@test sort([cancercell1A.genome.CN[1].tot, cancercell2A.genome.CN[1].tot]) == [1, 3]

cancercell1 = CINulator.copycell(cell)
cancercell2 = CINulator.copycell(cell)

cancercell1, cancercell2, Rmax, killcell1, killcell2, labelid =
        CINulator.newmutations(cancercell1, 
                                cancercell2, 
                                μ,
                                s, 
                                10.0,
                                0.0, 
                                CINulator.optimumfitness(),
                                maxCN = 6, 
                                minCN = 1,
                                labelcells = true,
                                labelid = 1)

@test (map(x -> x.tot, cancercell1.genome.CN) .+ map(x -> x.tot, cancercell2.genome.CN)) / 2 == fill(2.0, Nchr)


@time simresult = CINulator.simulate(b, d, 
                10^3, 
                Nchr,
                μ = μ, 
                s = s,
                fitnessfunc = CINulator.optimumfitness(),
                labelcells = true);

@time simresult = CINulator.simulate(b, d, 
                10^3, 
                Nchr,
                μ = μ, 
                s = s,
                fitnessfunc = CINulator.optimumfitness(),
                states = [(1,2), (2,3)],
                labelcells = true);


