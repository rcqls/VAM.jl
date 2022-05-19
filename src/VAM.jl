module VAM

export Model

include("compute.jl")
include("model.jl")
include("family.jl")
include("maintenance.jl")
include("stop_policy.jl")
include("simulate.jl")
include("mle.jl")

end
