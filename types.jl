module typs

@enum AmbStatus begin #Define the ambulace status types
    Idle
    Bussy
    Null
end

n_vertiports=2
n_eVTOL=4
eVTOL_speed=165

mutable struct Ambulance #Create the ambulance structure
    index::Int64
    region::Int64
    station::Int64
    class::Int64
    status::AmbStatus
    returnBase::Float64
    #Ambulance()=new(0,0,0,0,Idle,0.0)
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

export Ambulance, AmbStatus, Idle, Bussy, Null, Request2

end