
using CSV
using DataFrames
using DataStructures
using StatsBase
using Base.Math

@enum AmbStatus begin #Define the ambulace status types
    Idle
    Bussy
    Null
end

mutable struct Request2 # Create Request structure
    index::Int64
    region::Int64
    latitude::Float64
    longitude::Float64
    travtimesfrombases::Array{Float64}
    travttonearesth::Float64
    onsitetime::Float64
    traveltfromhtobases::Array{Float64}
    priority::Int64
    reqhospit::Int64
    inhardtoreachzone::Bool
    inqueue::Bool
    time_request::Float64
    response_time::Float64
    timetoHospital::Float64
    ServedbyAmbIndex::Int64
    finishedService::Float64
end


freq_calls=[[14,  10,  11,  11,  11,  13,  11],
            [9,  6,  7,  7,  7,  8,  7],
            [12,  9,  10,  10,  10,  12,  9],
            [7,  5,  6,  6,  6,  7,  6],
            [18,  13,  15,  15,  15,  17,  14],
            [4,  3,  4,  4,  4,  4,  3],
            [23,  16,  19,  19,  19,  21,  17],
            [53,  38,  44,  44,  44,  50,  41],
            [21,  15,  17,  17,  17,  19,  16],
            [11,  7,  9,  9,  9,  10,  8]]

freq_horaria_acc=[0.027,  0.049,  0.074,  0.096,  0.117,  0.1365, 
             0.165,  0.196,  0.243,  0.303,  0.365,  0.425,  0.481,  
             0.539,  0.59,  0.638,  0.691,  0.743,  0.797,  0.842,  0.892,  0.932,  0.97,  1]





             
#objetivo es seleccionar las requisiciones aleatoriamente con base en región, semana y hora


# at least one EVTOL per open vertiport n>=m


# crear un arreglo con los 

# 