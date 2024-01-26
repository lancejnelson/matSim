module dataset


include("/Users/legoses/OneDrive - BYU-Idaho/codes/Crystal.jl")
using .CrystalMethods: Crystal,getEnergy,writePOSCAR,CartesianToDirect!,DirectToCartesian
using Printf
using DelimitedFiles
using Statistics

struct DataSet
    crystals::Vector{Crystal}
end


function readStructuresIn(folder::String,file::String)
    cd(folder)
    file = open(file,"r")
    pos = readlines(file)

    data = Vector{Crystal}()
    for (idx,line) in enumerate(pos)

        if occursin("#--",line)
            nAtoms = sum([parse(Int64,x) for x in split(pos[idx + 6])])
            startpoint = idx + 1
            theend = idx + 7 + nAtoms
            thisCrystal = Crystal(pos[startpoint:theend],energyFP = parse(Float64,pos[theend + 2 ]))
            if !isnan(thisCrystal.energyFP)
                push!(data,thisCrystal)
            end
        end
    end
    return DataSet(data)
end


function writeStructuresIn(path::String, structures::DataSet)
    io = open(path, "w")
    for crystal in structures.crystals
        write(io,"#-----------------------\n")
        write(io,crystal.title * "\n")
        write(io,string(crystal.latpar) * "\n")
        for lv in eachcol(crystal.lVecs)
            write(io,join(lv," ") * "\n")
        end
        write(io,join(crystal.nType," ") * "\n")
        write(io,crystal.coordSys[1] * "\n")
        for typ in crystal.atomicBasis
            for atom in typ
                write(io,join(atom," ") * "\n")
            end
        end
        write(io,"#Energy:\n")
        write(io,string(crystal.energyFP) * "\n")
    end
    close(io)
    
end



function readVaspFolders(folder::String,file::String;poscar = "CONTCAR",outcar = "OUTCAR")

    for obj in readdir(folder,join = true)

        if isdir(obj)
            println(obj)
            poscarPresent = isfile(joinpath(obj,poscar))
            outcarPresent = isfile(joinpath(obj,outcar)) 
            if poscarPresent && outcarPresent
                crystal = Crystal(obj,poscar,energyFP = getEnergy(joinpath(obj,outcar)))
                writePOSCAR(crystal,joinpath(folder,file),"a")
                open(joinpath(folder,file), "a") do f
                    write(f,"Energy:","\n")
                    writedlm(f,crystal.energyFP)
                    write(f,"#--------------------","\n")
                end
                
                
            else
                @printf("Didn't find one of the necessary files. POSCAR: %s OUTCAR: %s Skipping this dir: %s\n",poscarPresent ? "true" : "false",outcarPresent ? "true" : "false", obj)
                #print(obj)
            end
        end
    end



end

end