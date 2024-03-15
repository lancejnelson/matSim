using Plots
using Statistics
using Printf
using Distributions

include("/Users/legoses/OneDrive - BYU-Idaho/codes/modules/LennardJones.jl")
using .LJ: readStructuresIn,totalEnergy,getTraining_Holdout_Sets,initializeMetrop,getSamples #,findFit Everything having to do with crystals... Reading from file.. Writing to file. Manipulating.. etc. 
dset = readStructuresIn("/Users/legoses/Downloads","structures.CuPt")
println(length(dset.crystals))
display(dset.crystals[4].atomicBasis)
## Standardize the data..
meanEnergy = mean([i.energyFP for i in dset.crystals])
stdEnergy = std([i.energyFP for i in dset.crystals])
offset = 3
for i in dset.crystals
    i.energyFP = (i.energyFP - meanEnergy)/stdEnergy-offset
end

## Pre-calculate the distances needed for LJ
cutoff = 20.6
for crystal in dset.crystals
    crystal.ljvals .= totalEnergy(crystal,cutoff)
end
# Split data set into training and holdout sets
trainingSet, holdoutSet = getTraining_Holdout_Sets(dset,150)


# Specify Metropolis settings.

nDraws = 10000
nBurnIn = 1000
candSig = [0.1 0.04
          0.08 0.04
          0.1 0.04]

muGuess = [0.5 1.2
          1.8 1.5
          1.3 1.4]

sigmaGuess = 2.5
candSig_sig = 0.01


#Run Metropolis-Hastings on training set.

μprior = [Gamma(3,2) Gamma(70,30)
         Gamma(3,2) Gamma(70,30)
         Gamma(3,2) Gamma(70,30) ]
σprior = Gamma(3,2)
logPost(data,μ,σ) = LJ.loglik(data,μ,σ) +  sum(logpdf.(μprior,μ)) + logpdf(σprior,σ)
proposal(μ,σ) = Gamma(μ^2/σ^2,σ^2/μ)

thisMetrop = initializeMetrop(trainingSet,nDraws,nBurnIn,candSig,candSig_sig,muGuess,sigmaGuess,proposal,logPost)
results = getSamples(thisMetrop)

# Look at the acceptance rates.
display(results.μ_accept)
println(results.σ_accept)

# Inspect the histograms of the LJ parameters
histogram(results.sigmadraws)
histogram(results.μ_draws[5*nBurnIn:end,2,1],bins = 100)
plot(results.μ_draws[5*nBurnIn:end,2,2],serietype= :scatter)
histogram([i.energyFP for i in dset.crystals])

# Predict on the holdout set to evaluate quality of model.
trueVals = [(i.energyFP + offset) * stdEnergy + meanEnergy for i in holdoutSet.crystals]
predictVals =  [mean([ (totalEnergy(j,results.mudraws[i,:,:]) + offset) *stdEnergy + meanEnergy  for i in 10000:100000]) for j in holdoutSet.crystals]
predictUnc =  [2 * std([(totalEnergy(j,results.mudraws[i,:,:]) + offset) *stdEnergy + meanEnergy  for i in 10000:100000]) for j in holdoutSet.crystals]

rmsError = sqrt(mean([(trueVals[i] - predictVals[i])^2 for i in 1:length(holdoutSet.crystals)]))
x = -7:0.05:-3

# Plot predicted vs true energies. Slope=1 line is a perfect fit.
myp = plot(predictVals,trueVals,seriestype = :scatter,xerror = predictUnc,ms = 2.5,ylabel = "True Energy (eVs/atom)", xlabel = "Predicted Energy (eVs/atom)",legend=false)
plot!(x,x,lw = 5)


tag = @sprintf("RMS Error: %5.2f",rmsError)
annotate!(-6,-3.2,tag)


savefig(myp,"agpt.png")

