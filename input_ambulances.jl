module inputamb

include("types.jl")

using .typs
using CSV
using DataFrames


function readAmbData(AmbFileName::String) #Read the ambulance file and return it into a dataframe
    dfambulances=DataFrame(CSV.File(AmbFileName))
    Ambulances=Array{Ambulance}(undef,size(dfambulances,1))
    for (i, row) in enumerate(eachrow(dfambulances))
        Ambulances[i]=Ambulance(dfambulances[i,:index],dfambulances[i,:region],dfambulances[i,:station],dfambulances[i,:class],Idle,0.0)
    end
    return Ambulances
end


function Inputambs(AmbFileName::String)
    ListAmb=readAmbData(AmbFileName)
    dfambulances2=DataFrame(CSV.File(AmbFileName))
    numReg=length(unique(dfambulances2.region))
    Num_estaciones=Vector(combine(groupby(dfambulances2, :region), :station => maximum => :station)[:, "station"])
    NumAmbulances=length(ListAmb)
    resources=Vector{Vector{Vector{Vector{Ambulance}}}}(undef, numReg)
    for i in 1:numReg
        resources[i] = Vector{Vector{Vector{Ambulance}}}(undef, Num_estaciones[i])
        for j in 1:Num_estaciones[i]
            resources[i][j] = [[],[]]
        end
    end
    for i in 1:NumAmbulances
    push!(resources[ListAmb[i].region][ListAmb[i].station][ListAmb[i].class],ListAmb[i])
    end
    return resources
end

function numRegiones(AmbFileName::String)
    dfambulances2=DataFrame(CSV.File(AmbFileName))
    numReg=length(unique(dfambulances2.region))
    return numReg
end

export Inputambs, numRegiones, Ambulance, AmbStatus, Idle, Bussy, Null, Request2

end