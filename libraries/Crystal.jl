

Crystal(folder::String,file::String; energyFP = 0, energyPred = 0) = fromPOSCAR(folder,file, energyFP = 0, energyPred = 0)
Crystal(list::Vector{String}; energyFP = 0, energyPred = 0) = fromPOSCAR(list, energyFP = energyFP, energyPred = energyPred)


function getPartitions(atoms::Vector{Vector{Float64}},nAtoms::Vector{Int64})
    if length(atoms) != sum(nAtoms)
        println("Number of atoms doesn't match")
        return
    end
    parts = prepend!(cumsum(nAtoms),0)
    return [atoms[parts[i]+1:parts[i+1]] for i=1:length(parts)-1]
end

function buildRandom(lPar:: Float64, lVecs::Matrix{Float64},nAtoms::Vector{Int64},cutoff::Float64)

    totalAtoms = sum(nAtoms)
    order = length(nAtoms)
    aType = hcat([ [n for i=1:nAtoms[n]]' for n=1:length(nAtoms)]...)'

    atoms = [zeros(Float64,3) for i=1:totalAtoms]  

    nTotal = 0
    counter = 0
    while nTotal < totalAtoms
        counter += 1
        if counter > 5000
            println("Having trouble putting that many atoms into this simulation cell.")
            println("So far I've only place $nTotal atoms." )
            return
        end
        newAtomCart = lPar * lVecs * rand(3)

        newOK = true
        for i=1:nTotal
            if norm(newAtomCart - atoms[i]) < cutoff
                newOK = false
            end
            if !newOK
                break
            end
        end
        if newOK
            nTotal += 1
            atoms[nTotal] .= newAtomCart
        end
    end
    atomicBasis = getPartitions(atoms,nAtoms)
    return Crystal("Random Locations",lPar,lVecs,nAtoms,aType,totalAtoms,["C"], atomicBasis,0.0,0.0,order,zeros(3,2))
end


function fromPOSCAR(folder::String,file::String;energyFP = 0, energyPred = 0)

    cd(folder)
    
    file = open(file, "r")
    pos = readlines(file)

    title = pos[1]
    lVecs = SMatrix{3,3}(reduce(hcat,[parse.(Float64, y) for y in [split(x) for x in pos[3:5]]])) # Read in the lattice vectors.
    if !isletter(lstrip(pos[6])[1])
        counter = 6
    else
        counter = 7
    end
    nBasis = parse.(Int,split(pos[counter]))
    coordSys = [pos[counter + 1]]
    latpar = parse(Float64,pos[2])
    aType = hcat([ [n for i=1:nBasis[n]]' for n=1:length(nBasis)]...)'
    println(pos[8:7 + sum(nBasis)])
    allBasis = [SVector{3,Float64}(parse.(Float64,split(x)[1:3])) for x in pos[(counter + 2):(counter + 1 +sum(nBasis))]] # Read all of the basis vectors
    allTypes = try
        [split(x)[end] for x in pos[(counter + 2):end]]
    catch e
        println("No atom types listed in the POSCAR")
        ["?" for x in pos[(counter + 2):end]]
    end
    atomicBasis = getPartitions(allBasis,nBasis)
#    if length(nBasis) == 3
#        partitionIndices = [[1,nBasis[1]],[nBasis[1] + 1,nBasis[1] + nBasis[2]],[nBasis[1] + nBasis[2]+1,sum(nBasis)]] 
#    elseif length(nBasis) == 2
#        partitionIndices = [[1,nBasis[1]],[nBasis[1] + 1,sum(nBasis)]]
#    else
#        partitionIndices = [[1,nBasis[1]]]
#    end
#    atomicBasis = [allBasis[n[1]:n[2]] for n in partitionIndices] # Partition the list according to atom type.
    order = length(nBasis)
    nAtoms = sum(nBasis)
    return Crystal(title, latpar,lVecs,nBasis,aType,nAtoms,coordSys,atomicBasis,energyFP,energyPred,order,zeros(3,2))  # Create new crystal object.


#    pos = read(path,String)
#    lines = split(pos,"\n")
#    title = String(lines[1])
#    latpar = parse(Float64,lines[2])
#    lVecs = hcat(([parse(Float64,x) for x in split(y) ] for y in lines[3:5])...)
#    if !isletter(lstrip(lines[6])[1])
#        counter = 6
#    else
#        counter = 7
#    end
#    nBasis = [parse(Int64,x) for x in split(lines[counter])]
#    nAtoms = sum(nBasis)
#    coordSys = [lines[counter + 1][1]]
#    anBasis = cumsum(hcat([0],transpose(nBasis))[:])
#    aBasis = [[parse.(Float64,split(y)[1:3]) for y in lines[(counter + 2 +anBasis[ba]):(counter + 2 + anBasis[ba+1]-1)]] for ba in 1:length(nBasis)]
#    return Crystal(title,latpar,lVecs,nBasis,nAtoms,coordSys,aBasis,energyFP/nAtoms,energyPred,zeros(3,2))
end


function fromPOSCAR(lines::Vector{String};energyFP = NaN,energyPred = NaN)
    title = lines[1]
    lVecs = SMatrix{3,3}(reduce(hcat,[parse.(Float64, y) for y in [split(x) for x in lines[3:5]]])) # Read in the lattice vectors.
    if !isletter(lstrip(lines[6])[1])
        counter = 6
    else
        counter = 7
    end
    nBasis = parse.(Int,split(lines[counter]))
    coordSys = [lines[counter + 1]]
    latpar = parse(Float64,lines[2])
    aType = hcat([ [n for i=1:nBasis[n]]' for n=1:length(nBasis)]...)'
    allBasis = [SVector{3,Float64}(parse.(Float64,split(x)[1:3])) for x in lines[(counter + 2):(counter + 1)+sum(nBasis)]] # Read all of the basis vectors
    allTypes = try
        [split(x)[end] for x in lines[(counter + 2):end]]
    catch e
        println("No atom types listed in the POSCAR")
        ["?" for x in pos[(counter + 2):end]]
    end
    if length(nBasis) == 3
        partitionIndices = [[1,nBasis[1]],[nBasis[1] + 1,nBasis[1] + nBasis[2]],[nBasis[1] + nBasis[2]+1,sum(nBasis)]] 
    elseif length(nBasis) == 2
        partitionIndices = [[1,nBasis[1]],[nBasis[1] + 1,sum(nBasis)]]
    else
        partitionIndices = [[1,nBasis[1]]]
    end
    atomicBasis = [allBasis[n[1]:n[2]] for n in partitionIndices] # Partition the list according to atom type.
    order = length(nBasis)
    nAtoms = sum(nBasis)
    return Crystal(title, latpar,lVecs,nBasis,aType,nAtoms,coordSys,atomicBasis,energyFP,energyPred,order,zeros(3,2))  # Create new crystal object.

end



function DirectToCartesian!(crystal::Crystal)
    if crystal.coordSys[1] == "D"
        println("Converting to cartesian")
        crystal.atomicBasis .= [[crystal.lVecs * i for i in j] for j in crystal.atomicBasis ]
        crystal.coordSys[1] = "C"
    else
        println("Already in Cartesian coordinates")
    end
end

function CartesianToDirect!(crystal::Crystal)
    if crystal.coordSys[1] == "C"
        println("Converting to direct")
        crystal.atomicBasis .= [[inv(crystal.lVecs) * i for i in j] for j in crystal.atomicBasis ]

        crystal.coordSys[1] = "D"
    else
        println("Already in Cartesian coordinates")
    end
end

function DirectToCartesian(lVecs::Matrix{Float64},atom:: Vector{Float64})
    return lVecs * atom
end

function CartesianToDirect(lVecs::Matrix{Float64},atom:: Vector{Float64})
    return inv(lVecs) * atom
end

function isEquivalentAtom(atomOne,atomTwo)  # Atoms are assumed to be in direct coordinates.

    #=If the atoms are separated by a lattice vector we can be sure that they are equivalent atoms (translationally).
       To check this, subtract the two vectors and see if the resulting vector is only composed of 
       integers. 
    =#
    return all(round.(Int64,atomOne - atomTwo) - (atomOne - atomTwo) == [0,0,0]) 


end


function writePOSCAR(crystal::Crystal,fName::String,openCode::String = "w")
    open(fName,openCode) do f
        write(f,crystal.title,'\n')
        writedlm(f,crystal.latpar',' ')
        writedlm(f,crystal.lVecs',' ')
        writedlm(f,crystal.nType',' ')
        write(f,crystal.coordSys[1],'\n')
        counter = 1
        for basis in crystal.atomicBasis
            for atom in basis
                writedlm(f,atom',' ')
               # seek(f,position(f) - 1)
            #    write(f,' ' * crystal.atomTypes[counter],'\n')
                counter +=1
            end
        end
    
    end
    

end

function getEnergy(filePath::String)

    open(filePath) do f
        lines = readlines(f)
        energyLine = findall([occursin("free  energy",x) for x in lines])
        return parse(Float64,split(lines[energyLine[1]])[end-1])
    end

end

