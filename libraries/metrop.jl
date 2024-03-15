

function initializeMetrop(data::DataSet,nDraws::Int, nBurnIn::Int, candSig:: Array{Float64,2},candSig_sig::Float64, μGuess::Array{Float64,2}, σGuess::Float64,proposal::Function,logpost::Function)
    nTotal = nBurnIn + nDraws
    nParams = size(candSig)
    μ_draws = zeros(nTotal,nParams...)
    μ_draws[1,:,:] .= μGuess
    σ_draws = zeros(nTotal)
    σ_draws[1] = σGuess
    μ_accept = zeros(nParams...)
    σ_accept = Float64[0]

    LJMetrop(data,nTotal,nBurnIn,μ_draws,σ_draws,candSig,candSig_sig,μGuess,σGuess,μ_accept,σ_accept, proposal,logpost)
end



function getSamples(metrop::LJMetrop)
    nParams = size(metrop.candSig_μ)
    drawsWithCand = zeros(nParams...)
    for i = 2:metrop.nDraws
        metrop.μ_draws[i,:,:] .= metrop.μ_draws[i-1,:,:]  # Set the next draw to be equal to the previous.  I
        metrop.σ_draws[i] = metrop.σ_draws[i-1]
        for j = 1: nParams[1], k = 1:nParams[2]
            cand = rand(metrop.proposal(metrop.μ_draws[i,j,k],metrop.candSig_μ[j,k]))
            if cand < 0.05
                continue
            end 

            drawsWithCand .= metrop.μ_draws[i,:,:]  # Need to assemble the vector of parameters with the candidate draw inserted at the right place.
            drawsWithCand[j,k] = cand
            numerator = metrop.logpost(metrop.data,drawsWithCand,metrop.σ_draws[i]) + log(pdf(metrop.proposal(cand,metrop.candSig_μ[j,k]),metrop.μ_draws[i,j,k]))
            denominator =metrop.logpost(metrop.data,metrop.μ_draws[i,:,:],metrop.σ_draws[i])  + log(pdf(metrop.proposal(metrop.μ_draws[i,j,k],metrop.candSig_μ[j,k]),cand))
            r = numerator - denominator
            unif = log(rand(Uniform(0,1)))  # Draw from a uniform.
            if r >= 0 || ((r < 0) & (unif < r))  # Accept?
                metrop.μ_draws[i,j,k] = cand   # Yes!
                metrop.μ_accept[j,k] += 1/metrop.nDraws
            end
        end


        ## Now get sigma draws...

        cand = rand(metrop.proposal(metrop.σ_draws[i],metrop.candSig_σ))
        if cand < 0.05
            continue
        end 
        numerator = metrop.logpost(metrop.data,metrop.μ_draws[i,:,:],cand) + log(pdf(metrop.proposal(cand,metrop.candSig_σ),metrop.σ_draws[i]))
        denominator =metrop.logpost(metrop.data,metrop.μ_draws[i,:,:],metrop.σ_draws[i])  + log(pdf(metrop.proposal(metrop.σ_draws[i],metrop.candSig_σ),cand))
        r = numerator - denominator
        unif = log(rand(Uniform(0,1)))  # Draw from a uniform.
        if r >= 0 || ((r < 0) & (unif < r))  # Accept?
            metrop.σ_draws[i] = cand   # Yes!
            metrop.σ_accept[1] += 1/metrop.nDraws
        end

    end     
    return metrop
end




#function getSamples(metrop::Metrop)
#    drawsWithCand = zeros(metrop.nParams)
#    for i = 1:metrop.nDraws
#        metrop.draws[i,:] .= metrop.draws[i-1,:]  # Set the next draw to be equal to the previous.  I
#        for j = 1: metrop.nParams
#            cand = rand(Normal(draws[i,j],metrop.candSig[j]))  # Get a candidate draw.  Draw from distribution with previous draw
#            drawsWithCand[:] .= metrop.draws[:]  # Need to assemble the vector of parameters with the candidate draw inserted at the right place.
#            drawsWithCand[j] = cand
#            r = metrop.logpost(data,drawsWithCand...) - logga(data,metrop.draws[i,:]...)  #Ratio between last draw with 
#            unif = log(rand(Uniform(0,1)))  # Draw from a uniform.
#            if r >= 0 || (r < 0 && unif < r)  # Accept?
#                metrop.draws[i,j] .= cand   # Yes!
#                metrop.accept[j] .+= 1/metrop.nDraws
#            end
#        end
#    end     
#    println("acceptance rates", metrop.accept * 100)
#    return metrop
#end


#end

#using Pkg
#Pkg.add("Plots")
#using Distributions
#using DelimitedFiles
#using Plots
#using SpecialFunctions
#
#
#cd("/Users/legoses/Downloads/")
#data = readdlm("projectile.csv",',',Float64) # Read in the data
#μ = mean(data[:,1])
#σ = std(data[:,1])
#data[:,1] .= (data[:,1] .- mean(data[:,1]))./std(data[:,1])
#data[:,2] .= (data[:,2] .- mean(data[:,2]))./std(data[:,2])
#draws = bigLoop(data)
#
#
#a = (draws[:,1] .- μ)/σ
#histogram(a)
#plot(draws[:,3])
#plot(data[:,1],data[:,2],seriestype = :scatter)
#x = 0:0.01:10
#x = (x .- mean(x))/std(x)
#y = [draws[i,1] .* x.^2 .+ draws[i,2] .* x .+ draws[i,3] for i in 20000:1000:nDraws]
#display(draws[:,1])
#plot(x,y,legend=false)
#histogram(draws[10000:end,4],bins = 200)
#histogram2d(draws[10000:end,1], draws[10000:end,2])
#
#maximum(draws[2000:,2])
#minimum(draws[:,2])
#xPred = 175
#yPred = [rand(Normal(draws[i,1] * xPred^2 + draws[i,2]*xPred + draws[i,3],sqrt(21))) for i in 2000:100000]
#histogram(yPred)
#
#
#a = 0:0.01:3
#y = [logga(data,x,4,15) for x in a]
#maximum(y)
#plot(a,exp.(y))
#
#
#