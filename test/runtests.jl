using CINulator
using Distributions
using Test
using Random
using Statistics

tests = ["sim", "optim"]

println("Running tests ...")

for t in tests
    fn = "test_$t.jl"
    println("* $fn ...")
    include(fn)
end
