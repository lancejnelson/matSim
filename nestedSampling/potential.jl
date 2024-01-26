#======================================================================
=                           Potential
=
=   Finds the total potential of a 3-dimensional system
======================================================================#

function U(S, L, e=1.0, s=1.0)
    # returns the total energy of a 3-dimensional configuration
    # using the Lennard-Jones potential
    total = 0

    # loop and add the energy from each unique pair of atoms
    N = length(S[:,1])
    for i in 1:N
        for j in 1:(i-1)
            x = S[j,1] - S[i,1]
            y = S[j,2] - S[i,2]
            z = S[j,3] - S[i,3]

            # implement periodic boundary conditions
            if x > L/2
                x -= L
            elseif x < -L/2
                x += L
            end
            if y > L/2
                y -= L
            elseif y < -L/2
                y += L
            end
            if z > L/2
                z -= L
            elseif z < -L/2
                z += L
            end

            r = x^2 + y^2 + z^2
                
            total += 4e*((s^2/r)^6 - (s^2/r)^3)
        end
    end

    return total
end
