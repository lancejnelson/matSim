using Plots
using Statistics
using StatsBase
using Printf
push!(LOAD_PATH, "/Users/legoses/OneDrive - BYU-Idaho/codes/")
include("LennardJones.jl")
using .LennardJones: readStructuresIn,LJ,totalEnergy,findFit,getTraining_Holdout_Sets # Everything having to do with crystals... Reading from file.. Writing to file. Manipulating.. etc. 
dset = readStructuresIn("/Users/legoses/Downloads","structures.AgPt")

## Standardize the data..
meanEnergy = mean([i.energyFP for i in dset.crystals])
stdEnergy = std([i.energyFP for i in dset.crystals])

for i in dset.crystals
    i.energyFP = (i.energyFP - meanEnergy)/stdEnergy
end

## Pre-calculate the distances needed for LJ
cutoff = 30.6
for crystal in dset.crystals
    crystal.ljvals .= LJ(crystal,cutoff)
end
# Specify Metropolis settings.
nDraws = 10000
candSig = [0.002 0.005
          0.002 0.07
          0.002 0.05]

muGuess = [0.5 1.2
          1.8 1.5
          1.3 1.4]

sigmaGuess = 2.5
candSig_sig = 0.01
μpriorParams = reshape(Float64[3 2; 3 2
                               3 2; 70 30
                               70 30; 70 30],3,2,2)
σpriorParams = Float64[3; 2]

# Split data set into training and holdout sets
trainingSet, holdoutSet = getTraining_Holdout_Sets(dset,800)

#Run Metropolis-Hastings on training set.
results = findFit(trainingSet,nDraws,candSig,candSig_sig,μpriorParams,σpriorParams,muGuess,sigmaGuess)

# Look at the acceptance rates.
display(results.μaccept)
println(results.σaccept)

# Inspect the histograms of the LJ parameters
histogram(results.sigmadraws)
histogram(results.mudraws[3000:10000,2,:],bins = 100)


# Predict on the holdout set to evaluate quality of model.
trueVals = [i.energyFP  *stdEnergy + meanEnergy for i in holdoutSet.crystals]
predictVals =  [mean([totalEnergy(j,results.mudraws[i,:,:])  * stdEnergy + meanEnergy  for i in 1000:10000]) for j in holdoutSet.crystals]
predictUnc =  [std([totalEnergy(j,results.mudraws[i,:,:])  * stdEnergy + meanEnergy  for i in 1000:10000]) for j in holdoutSet.crystals]

rmsError = sqrt(mean([(trueVals[i] - predictVals[i])^2 for i in 1:length(holdoutSet.crystals)]))
x = -6:0.05:-3

# Plot predicted vs true energies. Slope=1 line is a perfect fit.
plot(predictVals,trueVals,seriestype = :scatter,xerror = predictUnc)
plot!(x,x)
tag = @sprintf("RMS Error: %5.2f",rmsError)
annotate!(-4.0,-3.2,tag)
