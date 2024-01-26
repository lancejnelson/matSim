#======================================================================
=                           -Random walk-
=
=   Takes some sample on a random walk. At each step, it checks whether
=   or not the energy is lower than the energy of the previous configu-
=   ration, and only applies the step if it lowers the energy.
======================================================================#

using Random
include("potential.jl")

function walk(S, W, L)
    N = length(S[:,1])

    E  = U(S, L)           # gets the energy of the current sample

    # does W random walks
    for i in 1:W

        # loops through each particle
        for j in 1:N
            # loop until the energy is lower than E

            vec = (rand(3) .- 0.5)*0.1   # interval [-0.1, 0.1)
            S[j,:] .+= vec
            (S[j,:] .+= L) .%= L            # wrap-around boundary conditions

            if U(S, L) >= E         # compare new energy to previous
                S[j,:] .-= vec
            end
            (S[j,:] .+= L) .%= L            # wrap-around boundary conditions
        end
    end
    
    return S
end
