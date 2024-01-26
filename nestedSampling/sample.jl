#======================================================================
=                           -Create a sample-
=
=   This file create a new sample with some number of atoms and a given
=   density. The atoms have a uniformly random distribution within the
=   space.
======================================================================#
using Random

# function that makes a new configuration
function create(N, L)
    # initialize all the particles to random spots in the lattice
    S = zeros(N, 3)

    for i in 1:N
        S[i,:] = rand(3)*L
    end

    return S
end
