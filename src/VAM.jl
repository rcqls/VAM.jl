module VAM

using Random, DataFrames
export Model

abstract type AbstractModel end

include("compute.jl")
include("formula.jl")
include("family.jl")
include("maintenance.jl")
include("maintenance_policy.jl")
include("stop_policy.jl")
include("model.jl")
include("simulate.jl")
include("mle.jl")

end
