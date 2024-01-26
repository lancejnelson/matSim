#======================================================================
                                -Nested Sampling-
                        For the Lennard-Jones potential

    e   -   For the LJ potential--multiplier that specifies the 
            strength of the potential. Default value is 1.0.
    s   -   For the LJ potential--specifies where the potential is 0.
            Default value is 1.0.
    N   -   Number of atoms in a sample.
    K   -   Specifies number of samples kept track of.
    W   -   Length of random walk.
    Ve  -   Volume threshold.
======================================================================#

include("potential.jl")
include("randwalk.jl")
include("sample.jl")

# Program constants
s = 1.0
e = 1.0
N = 64                              # atoms per sample
L = 12s                             # length of the cube
K = 512                             # samples we keep track of
W = 1024                            # length of random walk
Ve = 1e-6                           # volume threshold

# program variables
samples = zeros(Float64, N, 3, K)
E = zeros(Float64, K)

# creates samples and find energies
for i in 1:K
    samples[:,:,i] = create(N, L)
    E[i] = U(samples[:, :, i], L)
end

# Keep track of energies and volumes
V = [K/(K + 1)]
Em = Float64[]

n = 1

# only stops when some percentage of the space is left
while V[end] > Ve
    Emax, loc = findmax(E)
    append!(Em, Emax)

    samples[:,:,loc] = walk(samples[:,:,loc], W, L)   # random walk
    E[loc] = U(samples[:,:,loc], L)

    global n += 1
    append!(V, (K/(K + 1))^n)
end

# here we want to create a new file, that has a unique name, and print all of the output to
# that file, not to the slurm file. Something like output1, output2, output3, ...

println("s:\t", s)
println("e:\t", e)
println("N:\t", N)
println("K:\t", K)
println("W:\t", W)
println("Ve:\t", Ve)
println()

# print energies and volumes
for i in 1:length(Em)
    println(Em[i])
end
for i in 1:length(V)
    println(V[i])
end
