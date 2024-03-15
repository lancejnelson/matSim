



mutable struct Crystal
    title::String
    latpar::Float64
    lVecs:: Matrix{Float64}
    nType:: Vector{Int64} #Number of each type of atom 
    aType:: Vector{Int64} # Integers representing the type for each basis atom
    nAtoms::Int64  # Total number of atoms in the unit cell
    coordSys::Vector{String} # 'D' for Direct or 'C' for cartesian
    atomicBasis:: Vector{Vector{Vector{Float64}}}  # List of all the atomic basis vector separated by type  
    energyFP:: Float64  # First principles energy
    energyPred:: Float64 # Model energy
    order::Int64 # binary, ternary, etc.
    ljvals:: Matrix{Float64}
end

struct DataSet
    crystals::Vector{Crystal}
end


struct NS
    L :: Float64  # Size of simulations cell (cube)
    K :: Int64  # Number of configurations in simulation
    W :: Int64  # Length of random walk.
    V_ϵ :: Float64 # Convergence threshold
    N :: Vector{Int64}  # Number of atoms per configuration
    configs :: Vector{Crystal}  # Keeps track of all of the configurations in the simulation
#    E :: Vector{Float64}  # Energy of each configuration    
end

struct LJMetrop
    data::DataSet
    nDraws:: Int
    nBurnIn:: Int
    μ_draws:: Array{Float64,3}
    σ_draws:: Vector{Float64}
    candSig_μ:: Array{Float64,2}
    candSig_σ:: Float64
    μ_guess:: Array{Float64,2}
    σ_guess:: Float64
    μ_accept:: Array{Float64,2}
    σ_accept:: Vector{Float64}
    proposal:: Function
    logpost:: Function
end