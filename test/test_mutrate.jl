Nchr = 2
mut = 1
μ = CINulator.Chrmutrate(Nchr, mutrates = fill(mut, Nchr), psplit = fill(0.0, Nchr))
mygenome = Genome(2)
mygenome.CN[2] = Chromosome(Arm(2,2), Arm(2,2))
s = Chrfitness(mygenome, alpha = 0.5)
b = log(2)
d = 0.1

cell = CINulator.cancercellCN(log(2), 0.0, log(2), 0.0, Float64[], [], CINulator.Genome(Nchr), Int64[], "temp")
cancercell1 = CINulator.copycell(cell)
cancercell2 = CINulator.copycell(cell)

cancercell1A, cancercell2A = CINulator.missegregation(cancercell1, cancercell2, 1, 0);
@test (map(x -> x.p.tot, cancercell1A.genome.CN) .+ map(x -> x.p.tot, cancercell2A.genome.CN)) / 2 == fill(2.0, Nchr)
@test (map(x -> x.q.tot, cancercell1A.genome.CN) .+ map(x -> x.q.tot, cancercell2A.genome.CN)) / 2 == fill(2.0, Nchr)
@test cancercell1A.genome.CN[1].p.tot != cancercell2A.genome.CN[1].p.tot
@test sum([cancercell1A.genome.CN[1].p.tot, cancercell1A.genome.CN[1].q.tot] .== [cancercell2A.genome.CN[1].p.tot, cancercell2A.genome.CN[1].q.tot]) == 0
@test sort([cancercell1A.genome.CN[1].p.tot, cancercell2A.genome.CN[1].p.tot]) == [1, 3]

#test for when only one arm is split
cancercell1 = CINulator.copycell(cell)
cancercell2 = CINulator.copycell(cell)
cancercell1A, cancercell2A = CINulator.missegregation(cancercell1, cancercell2, 1, 1);
@test sum([cancercell1A.genome.CN[1].p.tot, cancercell1A.genome.CN[1].q.tot] .== [cancercell2A.genome.CN[1].p.tot, cancercell2A.genome.CN[1].q.tot]) == 1

cancercell1 = CINulator.copycell(cell)
cancercell2 = CINulator.copycell(cell)

cancercell1A, cancercell2A, Rmax, killcell1, killcell2, labelid =
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
                                

μ = CINulator.Chrmutrate(Nchr, mutrates = fill(mut, Nchr), psplit = fill(1.0, Nchr))

cancercell1 = CINulator.copycell(cell)
cancercell2 = CINulator.copycell(cell)

cancercell1A, cancercell2A, Rmax, killcell1, killcell2, labelid =
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
println(cancercell1A.genome)
println(cancercell2A.genome)

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


