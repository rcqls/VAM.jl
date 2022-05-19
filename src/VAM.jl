module VAM

using Random, DataFrames
export Model

include("compute.jl")
include("family.jl")
include("maintenance.jl")
include("stop_policy.jl")
include("model.jl")
include("simulate.jl")
include("mle.jl")

end
