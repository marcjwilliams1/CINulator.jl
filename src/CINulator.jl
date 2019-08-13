module CINulator

using Distributions
using StatsBase
using Printf
using Statistics
using Random
using RecipesBase
using Reexport
using DataFrames
@reexport using Plots
import Plots: _cycle
#using Plots.PlotMeasures

import Base.show

export
    copynumberfrequency,
    simulate,
    sitefrequency,
    celldataframe,
    mergecelldataframe,
    optimumfitness,
    Chrfitness,
    Chrmutrate,
    samplecells,
    simulate_timeseries


### source files
include("simulate.jl")
include("plots.jl")
include("util.jl")

end # module
