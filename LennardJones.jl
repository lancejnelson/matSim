
module LennardJones
using LinearAlgebra
using StaticArrays
using Distributions


include("/Users/legoses/OneDrive - BYU-Idaho/codes/dataset.jl")
## When building a lj model, I'll need a crystal data type and I'll need to be able to read a data file.
using .dataset:Crystal,readStructuresIn,DataSet,CartesianToDirect!,DirectToCartesian
#include("/Users/legoses/OneDrive - BYU-Idaho/codes/metrop.jl")
#using .metrop


struct Metrop
    data::DataSet
    nDraws:: Int
    mudraws:: Array{Float64,3}
    sigmadraws:: Vector{Float64}
    candSig:: Matrix{Float64}
    candSig_sig:: Float64
    μpriorParams:: Array{Float64,3}
    σpriorParams:: Vector{Float64}
    MH::Bool
    μaccept:: Matrix{Float64}
    σaccept:: Vector{Float64}
end


function LJ(crystal::Crystal,cutoff::Float64)
    CartesianToDirect!(crystal)
    ljvals = zeros(3,2)  #Specific to binary material.  Needs generalized to n-ary case.
    
    loopBounds = convert.(Int64,cld.(cutoff ,[norm(x) for x in eachcol(crystal.latpar * crystal.lVecs)] ))
    
    # The outer two loops are to loop over different centering atoms.
    for (iCenter,centerAtomType) in enumerate(crystal.atomicBasis) 
        for centerAtom in centerAtomType
            centerAtomC = DirectToCartesian(crystal.latpar * crystal.lVecs,centerAtom)
            # These two loops are to loop over all possible neighbor atoms
            for (iNeighbor,aType) in enumerate(crystal.atomicBasis)  #Loop over the different atom types.

                for neighboratom in aType  # Loop over each atom per type.
                    # And these three inner loops are to find all of the periodic images of a neighboring atom.
                    for i = -loopBounds[1]:loopBounds[1], j = -loopBounds[2]:loopBounds[2], k= -loopBounds[3]:loopBounds[3]

                        newAtom = neighboratom + [i, j, k]
                        newCart = DirectToCartesian(crystal.latpar * crystal.lVecs,newAtom)
                        r = norm(newCart - centerAtomC) 
                        if r < cutoff && !isapprox(r,0)
                            #If the neighbor atom is inside the unit cell, then its going to be
                            # double counted at some point when we center on the other atom.  
                            # So we count it as half each time.
                            if all([i,j,k] .== 0 ) 
                                ljvals[iNeighbor + iCenter - 1,1] +=  1/2 * 1/r^6
                                ljvals[iNeighbor + iCenter - 1,2] +=  1/2 * 1/r^12
                            else
                                ljvals[iNeighbor + iCenter - 1,1] += 1/r^6
                                ljvals[iNeighbor + iCenter - 1,2] += 1/r^12
                            end
                        end
                    end
                end
            end
        end
    end
    return ljvals
end

function totalEnergy(crystal::Crystal,params)
    energy = sum(4 * params[1,:] .* params[2,:] .^6 .* crystal.ljvals[1,:] - 4 * params[1,:] .* params[2,:] .^12 .* crystal.ljvals[2,:] )
    return energy
end


function logpostnorm_norm(data::DataSet,ljparams,σ,μ_ϵ,σ_ϵ)  # Normal-Normal Relationship
    n = length(data)
    #- (c-mu_c)^2/2/sigma_c^2
    thesum =  - 1/(2 * σ) * sum([(i.energyFP - totalEnergy(i,ljparams))^2   for i in data.crystals])- (ϵ-mu_ϵ)^2/2/sigma_ϵ^2
    return thesum

end

function findFit(data::DataSet,nDraws::Int, candSig:: Matrix{Float64},candSig_sig::Float64, μpriorParams:: Array{Float64,3},σpriorParams:: Vector{Float64}, muGuess::Matrix{Float64},sigmaGuess::Float64)
    nParams = size(candSig)
    #fpenergies = [i.energyFP for i in data.crystals]
    mudraws = zeros(nDraws,nParams...)
    mudraws[1,:,:] .= muGuess
    sigmadraws = zeros(nDraws)
    sigmadraws[1] = sigmaGuess
    μaccept = zeros(nParams...)
    σaccept = 0
    mymetrop = Metrop(data,nDraws,mudraws,sigmadraws,candSig,candSig_sig,μpriorParams,σpriorParams,false,μaccept,[σaccept])
    
    return getSamples(mymetrop)


end

function getSamples(metrop::Metrop)
    nParams = size(metrop.candSig)
    drawsWithCand = zeros(MMatrix{nParams...})
    for i = 2:metrop.nDraws
        metrop.mudraws[i,:,:] .= metrop.mudraws[i-1,:,:]  # Set the next draw to be equal to the previous.  I
        metrop.sigmadraws[i] = metrop.sigmadraws[i-1]
        for j = 1: nParams[1], k = 1: nParams[2]
            # These two lines are to get the gamma centered at the right location with the right variance.
            α = metrop.mudraws[i,j,k]^2/metrop.candSig[j,k]
            θ = metrop.candSig[j,k]/metrop.mudraws[i,j,k]
            cand = rand(Gamma(α,θ))  # Get a candidate draw.  Draw from distribution with previous draw
            if cand < 0.05
                continue
            end 
            α_2 = cand^2/metrop.candSig[j,k]
            θ_2 = metrop.candSig[j,k]/cand

            drawsWithCand[:,:] .= metrop.mudraws[i,:,:]  # Need to assemble the vector of parameters with the candidate draw inserted at the right place.
            drawsWithCand[j,k] = cand
            r = logpostμ(metrop.data,metrop.sigmadraws[i], drawsWithCand,(j,k),metrop.μpriorParams[j,k,:]) + log(pdf(Gamma(α_2,θ_2 ),metrop.mudraws[i,j,k])) - logpostμ(metrop.data,metrop.sigmadraws[i], metrop.mudraws[i,:,:],(j,k),metrop.μpriorParams[j,k,:]) - log(pdf(Gamma(α,θ),cand)) #Ratio between last draw with 
            unif = log(rand(Uniform(0,1)))  # Draw from a uniform.
            if r >= 0 || (r < 0 && unif < r)  # Accept?
                metrop.mudraws[i,j,k] = cand   # Yes!
                metrop.μaccept[j,k] += 1/metrop.nDraws
            end
        end

        ## Now get sigma draws...

        α = metrop.sigmadraws[i]^2/metrop.candSig_sig
        θ = metrop.candSig_sig/metrop.sigmadraws[i]
        cand = rand(Gamma(α,θ))  # Get a candidate draw.  Draw from distribution with previous draw
        if cand < 0.05
            continue
        end 
        α_2 = cand^2/metrop.candSig_sig
        θ_2 = metrop.candSig_sig/cand
        r = logpostσ(metrop.data,cand, metrop.mudraws[i,:,:],metrop.σpriorParams) + log(pdf(Gamma(α_2,θ_2 ),metrop.sigmadraws[i])) - logpostσ(metrop.data,metrop.sigmadraws[i], metrop.mudraws[i,:,:],metrop.σpriorParams) - log(pdf(Gamma(α,θ),cand)) #Ratio between last draw with 
        unif = log(rand(Uniform(0,1)))  # Draw from a uniform.
        if r >= 0 || (r < 0 && unif < r)  # Accept?
            metrop.sigmadraws[i] = cand   # Yes!
            metrop.σaccept[1] += 1/metrop.nDraws
        end
    end     
#    println("acceptance rates", metrop.μaccept * 100)
    return metrop
end

# Normal-Gamma on all parameters in LJ
function logpostμ(data::DataSet,σ,μ,paramOfInterest,priorParameters)
    n = length(data.crystals)
    #- (c-mu_c)^2/2/sigma_c^2
    α = priorParameters[1]
    β = priorParameters[2]
    x = μ[paramOfInterest...]
#    println("------")
#    println(paramOfInterest)
#    println(x)
#    println(log(x))
 #   println(σ)
#    println([i.energyFP for i in data.crystals])
#    println([totalEnergy(i,μ) for i in data.crystals])
    thesum =  - 1/(2 * σ) * sum([(i.energyFP - totalEnergy(i,μ))^2   for i in data.crystals]) + (α - 1) * log(x) - β * x
    return thesum


    return lognorm(data,σ, totalEnergy(μ)) + lognorm(μ[paramOfInterest],priorParameters...)


end

# Normal-Gamma on Sigma
function logpostσ(data::DataSet,σ,μ,priorParameters)
    n = length(data.crystals)
    α = priorParameters[1]
    β = priorParameters[2]
    #- (c-mu_c)^2/2/sigma_c^2
    thesum =  -n/2 * log(σ) - 1/(2 * σ) * sum([(i.energyFP - totalEnergy(i,μ))^2   for i in data.crystals])+ (α - 1) * log(σ) - β * σ
    return thesum

    return lognorm(data,σ, totalEnergy(μ)) + lognorm(μ[paramOfInterest],priorParameters...)


end

end
