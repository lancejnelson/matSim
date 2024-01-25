module CrystalMethods

using Printf
using DelimitedFiles
using StaticArrays

export Crystal,readStructuresIn,writeStructuresIn,CartesianToDirect!,DirectToCartesian,readVaspFolders,DataSet

struct Crystal
    title::String
    latpar::Float64
    lVecs:: SMatrix{3,3,Float64,9}
    nType:: Vector{Int64} #Number of each type of atom 
    aType:: Vector{Int64} # Integers representing the type for each basis atom
    nAtoms::Int64  # Total number of atoms in the unit cell
    coordSys::Vector{String} # 'D' for Direct or 'C' for cartesian
    atomicBasis:: Vector{Vector{SVector{3,Float64}}}  # List of all the atomic basis vector separated by type  
    energyFP:: Float64  # First principles energy
    energyPred:: Float64 # Model energy
    order::Int64 # binary, ternary, etc.
    ljvals:: Matrix{Float64}
end

#struct DataSet
#    crystals::Vector{Crystal}
#end



Crystal(folder::String,file::String; energyFP = 0, energyPred = 0) = fromPOSCAR(folder,file, energyFP = 0, energyPred = 0)
Crystal(list::Vector{String}; energyFP = 0, energyPred = 0) = fromPOSCAR(list, energyFP = energyFP, energyPred = energyPred)

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

#function readStructuresIn(folder::String,file::String)
#    cd(folder)
#    file = open(file,"r")
#    pos = readlines(file)
#
#    data = Vector{Crystal}()
#    for (idx,line) in enumerate(pos)
#
#        if occursin("#--",line)
#            nAtoms = sum([parse(Int64,x) for x in split(pos[idx + 6])])
#            startpoint = idx + 1
#            theend = idx + 7 + nAtoms
#            thisCrystal = Crystal(pos[startpoint:theend],energyFP = parse(Float64,pos[theend + 2 ]))
#            push!(data,thisCrystal)
#        end
#    end
#    return DataSet(data)
#end
function DirectToCartesian!(crystal::Crystal)
    if crystal.coordSys[1] == "D"
        println("Converting to cartesian")
        crystal.atomicBasis .= [[crystal.lVecs * i for i in j] for j in crystal.atomicBasis ]

#        for (i,aType) in enumerate(crystal.basis)
#            for (j,atom) in enumerate(aType)
#                crystal.basis[i][j] =  crystal.lVecs * atom
#            end
#
#        end
        crystal.coordSys[1] = "C"
    else
        println("Already in Cartesian coordinates")
    end
end

function CartesianToDirect!(crystal::Crystal)
    if crystal.coordSys[1] == "C"
        println("Converting to direct")
        crystal.atomicBasis .= [[inv(crystal.lVecs) * i for i in j] for j in crystal.atomicBasis ]
#        for (i,aType) in enumerate(crystal.basis)
#            for (j,atom) in enumerate(aType)
#                crystal.basis[i][j] =  inv(crystal.lVecs) * atom
#            end
#
#        end
        crystal.coordSys[1] = "D"
    else
        println("Already in Cartesian coordinates")
    end
end

function DirectToCartesian(lVecs::SMatrix{3,3,Float64,9},atom:: SVector{3,Float64})
    return lVecs * atom
end

function CartesianToDirect(lVecs::SMatrix{3,3,Float64,9},atom:: SVector{3,Float64})
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
    println(openCode)
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

#function writeStructuresIn(path::String, structures::DataSet)
#    io = open(path, "w")
#    for crystal in structures.crystals
#        write(io,"#-----------------------\n")
#        write(io,crystal.title * "\n")
#        write(io,string(crystal.latpar) * "\n")
#        for lv in eachcol(crystal.lVecs)
#            write(io,join(lv," ") * "\n")
#        end
#        write(io,join(crystal.nType," ") * "\n")
#        write(io,crystal.coordSys[1] * "\n")
#        for typ in crystal.atomicBasis
#            for atom in typ
#                write(io,join(atom," ") * "\n")
#            end
#        end
#        write(io,"#Energy:\n")
#        write(io,string(crystal.energyFP) * "\n")
#    end
#    close(io)
#    
#end

#function readVaspFolders(folder::String,file::String;poscar = "CONTCAR",outcar = "OUTCAR")
#
#    for obj in readdir(folder,join = true)
#
#        if isdir(obj)
#            println(obj)
#            poscarPresent = isfile(joinpath(obj,poscar))
#            outcarPresent = isfile(joinpath(obj,outcar)) 
#            if poscarPresent && outcarPresent
#                crystal = Crystal(obj,poscar,energyFP = getEnergy(joinpath(obj,outcar)))
#                writePOSCAR(crystal,joinpath(folder,file),"a")
#                open(joinpath(folder,file), "a") do f
#                    write(f,"Energy:","\n")
#                    writedlm(f,crystal.energyFP)
#                    write(f,"#--------------------","\n")
#                end
#                
#                
#            else
#                @printf("Didn't find one of the necessary files. POSCAR: %s OUTCAR: %s Skipping this dir: %s\n",poscarPresent ? "true" : "false",outcarPresent ? "true" : "false", obj)
#                #print(obj)
#            end
#        end
#    end
#
#
#end

end