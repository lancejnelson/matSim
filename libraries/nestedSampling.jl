function initializeSimulation(lPar,lVecs,minSep,nAtoms,nConfigs,nWalks,V_ϵ)
    configs = [buildRandom(lPar,lVecs,nAtoms,minSep) for i in 1:nConfigs]
    return NS(lPar,nConfigs,nWalks,V_ϵ,nAtoms,configs)
end

function randomWalk(config::Crystal,nWalk::Int64,LJcutoff::Float64,LJparams::Array{Float64,2})
    # Loop over the number of random walks to take.
    for iWalk in 1:nWalk
        #Loop over all of the atoms in the simulation.
        for iType in 1:config.order, iAtom in 1:config.nType[iType]
            oldEnergy = sum(totalEnergy(config,LJcutoff,params = LJparams))
            newEnergy = oldEnergy + 100
            # Loop until you find a displacement that is downhill
            while newEnergy > oldEnergy
                randDisplacement = (rand(3) .- 0.5)*0.1 # Get a random displacement vector between -0.1 and 0.1
                config.atomicBasis[iType][iAtom] .+= randDisplacement  # Use the displacement to move this atom.
                newEnergy = sum(totalEnergy(config,LJcutoff,params = LJparams))
            end

        end
    end
end
