module CINulator

using Distributions
using StatsBase
using Printf
using Statistics
using Random
using RecipesBase
using Reexport
@reexport using Plots
import Plots: _cycle
#using Plots.PlotMeasures

import Base.show

export
    copynumberfrequency,
    simulate,
    sitefrequency


### source files
include("simulate.jl")
include("plots.jl")

end # module
