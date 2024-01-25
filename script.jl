using Plots
push!(LOAD_PATH, "/Users/legoses/OneDrive - BYU-Idaho/codes/")
include("LennardJones.jl")
using .LennardJones: readStructuresIn,LJ,totalEnergy,findFit # Everything having to do with crystals... Reading from file.. Writing to file. Manipulating.. etc. 
dset = readStructuresIn("/Users/legoses/Downloads","structures.AgPt")
cutoff = 7.6
for crystal in dset.crystals
    crystal.ljvals .= LJ(crystal,cutoff)
end

# Specify Metropolis settings.
nDraws = 10000
candSig = [1.5 0.1
          0.2 0.01
          2.9 2.9]

muGuess = [0.5 1.2
          1.8 1.5
          1.3 1.4]

sigmaGuess = 2.5
candSig_sig = 0.8
μpriorParams = reshape(Float64[3 2; 3 2
                               3 2; 3 2
                               3 2; 3 2],3,2,2)
σpriorParams = Float64[3; 2]

results = findFit(dset,nDraws,candSig,candSig_sig,μpriorParams,σpriorParams,muGuess,sigmaGuess)

display(results.μaccept)
println(results.σaccept)
histogram(results.sigmadraws)
histogram(results.mudraws[1000:10000,3,2],bins = 50)
