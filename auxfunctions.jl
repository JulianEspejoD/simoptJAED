module auxfunct

include("types.jl")

using .typs

function createqueues(numReg)
    queuesRegions2=Vector{Vector{Vector{Main.inputreq.typs.Request2}}}(undef, numReg)
    for i in 1:numReg
        queuesRegions2[i] = Vector{Vector{Main.inputreq.typs.Request2}}(undef, 3)
        for j in 1:3
            queuesRegions2[i][j] = []
        end
    end
    return queuesRegions2
end 

function randomdecisions(numReg, n_vertiports, n_evtols, Num_estaciones, dfLocation)
    auxarrayregions = Vector{Vector{Int64}}(undef, numReg)
    for i in 1:numReg
        auxarrayregions[i] = Int64[]  
        for j in 1:Num_estaciones[i]
            push!(auxarrayregions[i], j)
        end
    end
    Localizationvertiports= Array{Tuple{Int, Int}}(undef, n_vertiports)
    for i in 1:n_vertiports
        indices_nonEmptyreg = [i for (i, arr) in enumerate(auxarrayregions) if !isempty(arr)]
        randregion=indices_nonEmptyreg[rand(1:length(indices_nonEmptyreg))]
        randstation=rand(1:length(auxarrayregions[randregion]))
        Localizationvertiports[i]=(randregion,auxarrayregions[randregion][randstation])
        deleteat!(auxarrayregions[randregion],randstation)
    end
    List_evtols=Array{eVTOL}(undef,n_evtols)
    for i in 1:n_vertiports
    List_evtols[i]=eVTOL(1000+i,Localizationvertiports[i][1],Localizationvertiports[i][2],0,Idle,0.0,0.0,0.0,0.0)
    end
    for i in (n_vertiports+1):n_evtols
        randvertiport=rand(1:n_vertiports)
        List_evtols[i]=eVTOL(1000+i,Localizationvertiports[randvertiport][1],Localizationvertiports[randvertiport][2],0,Idle,0.0,0.0,0.0,0.0)
    end
    for element in List_evtols
        row=findfirst((dfLocation.region .== element.region) .& (dfLocation.station .== element.station))
        element.relative_pos=row
        element.latitude=dfLocation.latitude[row]
        element.longitude=dfLocation.longitude[row]
    end
    return List_evtols
end


function haversine(lat1, lon1, lat2, lon2)
    # distance between latitudes and longitudes
    dLat = (lat2 - lat1) * π / 180.0
    dLon = (lon2 - lon1) * π / 180.0
    # convert to radians
    lat1 = lat1 * π / 180.0
    lat2 = lat2 * π / 180.0
    # apply formulae
    a = sin(dLat / 2)^2 + sin(dLon / 2)^2 * cos(lat1) * cos(lat2)
    rad = 6371
    c = 2 * asin(sqrt(a))
    return rad * c
end




export createqueues, randomdecisions, haversine

end



#FUNCION PARA GUARDAR ESTADÍSTICAS DEL SISTEMA
