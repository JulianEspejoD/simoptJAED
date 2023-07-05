
using DiscreteEvents, Random, Distributions
using CSV
using DataFrames
using DataStructures
using StatsBase

include("input_ambulances.jl")
include("input_requests.jl")
include("auxfunctions.jl")
using .inputamb
using .inputreq
using .auxfunct

#INPUT DE AMBULANCIAS TERRESTRES
AmbFileName="inputAmbulances2.csv"
dfambulances2=DataFrame(CSV.File(AmbFileName))
numReg=length(unique(dfambulances2.region))
Num_estaciones=Vector(combine(groupby(dfambulances2, :region), :station => maximum => :station)[:, "station"])
NumAmbulances=size(dfambulances2)[1]
resourcess=inputamb.Inputambs2(dfambulances2, numReg, Num_estaciones, NumAmbulances)


#FUNCION PARA CREAR LAS COLAS DE REQUISICIONES SE PONE ACÁ PARA EVITAR CONFLICTOS
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
queuesRegions2=createqueues(numReg)


#INPUT DE EVTOL AMBULANCIAS 
n_vertiports=3
n_evtols=5
max_autnonomy=60
eVTOL_speed=165
LocationsFileName="input_location_samu.csv"
dfLocation=DataFrame(CSV.File(LocationsFileName))
ListeVTOL=randomdecisions(numReg, n_vertiports, n_evtols, Num_estaciones, dfLocation)
regs_eVTOL = unique([ev.region for ev in ListeVTOL])

#INPUT DE REQUISICIONES ENTRADA DE LLAMADAS

freq_calls=[[14,  10,  11,  11,  11,  13,  11],  #Frecuencia de llamadas por día de la semana
            [9,  6,  7,  7,  7,  8,  7],
            [12,  9,  10,  10,  10,  12,  9],
            [7,  5,  6,  6,  6,  7,  6],
            [18,  13,  15,  15,  15,  17,  14],
            [4,  3,  4,  4,  4,  4,  3],
            [23,  16,  19,  19,  19,  21,  17],
            [53,  38,  44,  44,  44,  50,  41],
            [21,  15,  17,  17,  17,  19,  16],
            [11,  7,  9,  9,  9,  10,  8]]

freq_horaria_acc=[0.027,  0.049,  0.074,  0.096,  0.117,  0.1365, #Probabilidad de llamadas durante el día Cumulative distribution function
             0.165,  0.196,  0.243,  0.303,  0.365,  0.425,  0.481,  
             0.539,  0.59,  0.638,  0.691,  0.743,  0.797,  0.842,  0.892,  0.932,  0.97,  1]
weeks_sim=8 #Cuantas semanas vamos a simular
ReqFileName="test_1_req.csv"

ListReq=readReqData2(ReqFileName)
arrivalsF=generatearrivals(ListReq, freq_calls, freq_horaria_acc, weeks_sim, numReg)
dim1array = reduce(vcat, arrivalsF)
sorted_arrivals = sort(dim1array, by = x -> x.time_request)

printing = true
global num_call=1

#ACA VAMOS A CREAR UN ARREGLO PARA COLECCIONAR ESTADÍSTICAS DEL SISTEMA

ReqServed=Array{Main.inputreq.typs.Request2,1}(undef, 0)
function statisticsReq(ServedReq::Main.inputreq.typs.Request2)
    push!(ReqServed,ServedReq)
end


#CUAL ES LA AMBUALNCIA AÉREA MÁS CERCANA


#cuantos vertipuertos hay en esa región
#ordenarlos de más cercano a más lejano
# filtrar los que estén en el rango
# mirar si están disponibles, si no hay ninguno se va con normal


function selectNearesteVTOL(ListeVTOL,resourcess,order_stations,reqArrived)
    evtol_in_reg = filter(x -> x.region == reqArrived.region, ListeVTOL)
    for element in evtol_in_reg
        element.auxdisthav = haversine(reqArrived.latitude, reqArrived.longitude, element.latitude, element.longitude)
    end
    evtol_in_reg_in_range=filter(x -> x.auxdisthav <= max_autnonomy, evtol_in_reg)
    if isempty(evtol_in_reg_in_range)
        return selectNearestAmbAv(resourcess,order_stations,reqArrived)   
    else
        distances = [element.auxdisthav for element in evtol_in_reg_in_range]
        perm = sortperm(distances)
        sorted_evtol_in_reg_in_range = evtol_in_reg_in_range[perm]
        for i in sorted_evtol_in_reg_in_range
            if i.status==Main.auxfunct.typs.Idle
                return i
            end     
        end
        return selectNearestAmbAv(resourcess,order_stations,reqArrived)  
    end
end

function selectNearestAmbAv(resourcess,order_stations,reqArrived)
    for k in order_stations #ACA TOCA ITERAR PRIMERO EN LA ESTACION estacion es 2
        for l in 1:length(resourcess[reqArrived.region][k][1])
            if (resourcess[reqArrived.region][k][1][l].status==Main.inputamb.typs.Idle)
                return resourcess[reqArrived.region][k][1][l]
             end
        end
        for l in 1:length(resourcess[reqArrived.region][k][2])
            if (resourcess[reqArrived.region][k][2][l].status==Main.inputamb.typs.Idle)
                return resourcess[reqArrived.region][k][2][l]
            end
        end
    end
    return nothing
end

function requestarrival(num_call)
    reqArrived=sorted_arrivals[num_call] # selecciona request aleatorio dentro de las llamadas
    global num_call+=1
    event!(𝐶,fun(requestarrival, num_call), at, sorted_arrivals[num_call].time_request) # programar el siguiente request
    order_stations=sortperm(reqArrived.travtimesfrombases)
    if reqArrived.inhardtoreachzone==true && reqArrived.region in regs_eVTOL
        amb_selected=selectNearesteVTOL(ListeVTOL,resourcess,order_stations,reqArrived)
    else
        amb_selected=selectNearestAmbAv(resourcess,order_stations,reqArrived)  #SELECCIONAMOS THE NEAREST AMBULANCE AVAILABLE
    end
    printing ? println(tau(), ": Request $(reqArrived.index) has arrived and the ambulance selected is $(amb_selected)") : nothing
    if amb_selected===nothing
        reqArrived.inqueue=true
        push!(queuesRegions2[reqArrived.region][reqArrived.priority],reqArrived)  #Coloca en la cola el request
    else
        serveEmergence(amb_selected,reqArrived)
    end
end


function serveEmergence(SelectedAmb,Request)# A ESTA FUNCIÓN LE MANDAMOS LA AMBULANCIA SELECCIONADA Y EL TRABAJO
    if isa(SelectedAmb, Main.auxfunct.typs.eVTOL)
        SelectedAmb.status=Main.auxfunct.typs.Bussy
        tt_evtol=(haversine(Request.latitude, Request.longitude, SelectedAmb.latitude, SelectedAmb.longitude)/eVTOL_speed)*60
        if (Request.inqueue==false)
            Request.response_time=Request.time_request+6+tt_evtol
        else
            Request.response_time=SelectedAmb.returnBase+6+tt_evtol    
        end
        Request.ServedbyAmbIndex=SelectedAmb.index
        printing ? println("Request $(Request.index) has been served by EVTOL: $(SelectedAmb.index)" ) : nothing
        return_to_base=Request.response_time+Request.onsitetime+tt_evtol
        Request.finishedService=return_to_base
        SelectedAmb.returnBase=return_to_base
        statisticsReq(Request)
        event!(𝐶, fun(followingReq,SelectedAmb), at, return_to_base)
    else    
        SelectedAmb.status=Main.inputamb.typs.Bussy
        Request.ServedbyAmbIndex=SelectedAmb.index
        printing ? println("Request $(Request.index) has been served by Ambulance: $(SelectedAmb.index)" ) : nothing
        if (Request.inqueue==false)
            Request.response_time=4+Request.travtimesfrombases[SelectedAmb.station]+Request.time_request #calculo de tiempo de respuesta aca llegamos al sitio
        else
            Request.response_time=4+Request.travtimesfrombases[SelectedAmb.station]+SelectedAmb.returnBase
        end
        #statisticsReq(Request)
        #printing ? println("Request $(Request.index) has been served at: $(Request.response_time) by Ambulance: $(SelectedAmb.index)" ) : nothing
        #aca revisamos si tenemos que llevarlo o no al hospital, si es tipo 1 puede si es 2 no
        if  (SelectedAmb.class==2) 
            return_to_base=Request.response_time+Request.onsitetime+Request.traveltfromhtobases[SelectedAmb.station] # TIEMPO DE RETORNO A BASE (response time mas + onsite time + devolucion a base)
            SelectedAmb.returnBase=return_to_base
            Request.finishedService=return_to_base
            statisticsReq(Request)
            event!(𝐶, fun(followingReq,SelectedAmb), at, return_to_base) # ACA PROGRAMAMOS EL EVENTO CUANDO VUELVA A BASE
            #printing ? println(return_to_base, ": Ambulance $(SelectedAmb.index) has returned to base") : nothing
            return nothing
        else
            if(Request.reqhospit==1)
                Request.timetoHospital=Request.response_time+Request.onsitetime+Request.travttonearesth
                return_to_base=Request.response_time+Request.onsitetime+Request.travttonearesth+ Request.traveltfromhtobases[SelectedAmb.station]# response time+ onsite time + travel to hospital time + return from hospital to base
                SelectedAmb.returnBase=return_to_base
                Request.finishedService=return_to_base
                statisticsReq(Request)
                event!(𝐶, fun(followingReq,SelectedAmb), at, return_to_base) # ACA PROGRAMAMOS EL EVENTO CUANDO VUELVA A BASE
                #printing ? println(return_to_base, ": Ambulance $(SelectedAmb.index) has returned to base") : nothing
                return nothing
            else
                return_to_base=Request.response_time+Request.onsitetime+Request.traveltfromhtobases[SelectedAmb.station] # TIEMPO DE RETORNO A BASE (response time mas + onsite time + devolucion a base)
                SelectedAmb.returnBase=return_to_base
                Request.finishedService=return_to_base
                statisticsReq(Request)
                event!(𝐶, fun(followingReq,SelectedAmb), at, return_to_base) # ACA PROGRAMAMOS EL EVENTO CUANDO VUELVA A BASE
                #printing ? println(return_to_base, ": Ambulance $(SelectedAmb.index) has returned to base") : nothing
                return nothing
            end
        end
    end
end


function priority(region, queuesRegions2)
    for queue in queuesRegions2[region]
        if !isempty(queue)
            return queue[1]
        end
    end
end

function priorityeVTOL(SelectedAmb,queuesRegions2)
    for queue in queuesRegions2[SelectedAmb.region]
        for req in queue
            dist=haversine(SelectedAmb.latitude, SelectedAmb.longitude, req.latitude, req.longitude)
            if dist <= max_autnonomy && req.inhardtoreachzone==true
                return req
            end
        end
    end
    return nothing
end


function followingReq(SelectedAmb)# #aca establece política de prioridad!!! EN ESTA FUNCION CAMBIAMOS EL STATUS DE LA AMBULANCIA 
    if isa(SelectedAmb, Main.auxfunct.typs.eVTOL)
        if (reduce(+, length.(queuesRegions2[SelectedAmb.region])) > 0)
            followingRequest=priorityeVTOL(SelectedAmb,queuesRegions2)
            if !isnothing(followingRequest)
                popfirst!(queuesRegions2[SelectedAmb.region][followingRequest.priority])
                event!(𝐶, fun(serveEmergence,SelectedAmb,followingRequest),at ,SelectedAmb.returnBase)
                printing ? println("EVTOL $(SelectedAmb.index) has returned to base at $(SelectedAmb.returnBase) and is going to: $(followingRequest.index)" ) : nothing
            else
                SelectedAmb.status=Main.auxfunct.typs.Idle     
            end                    
        else
            SelectedAmb.status=Main.auxfunct.typs.Idle
            printing ? println("EVTOL $(SelectedAmb.index) has returned to base at $(SelectedAmb.returnBase) and is going IDLE" ) : nothing
        end
    else
        if (reduce(+, length.(queuesRegions2[SelectedAmb.region])) > 0)
            followingRequest=priority(SelectedAmb.region,queuesRegions2)
            popfirst!(queuesRegions2[SelectedAmb.region][followingRequest.priority])
            event!(𝐶, fun(serveEmergence,SelectedAmb,followingRequest),at ,SelectedAmb.returnBase)
            printing ? println("Ambulance $(SelectedAmb.index) has returned to base at $(SelectedAmb.returnBase) and is going to: $(followingRequest.index)" ) : nothing
        else
            SelectedAmb.status=Main.inputamb.typs.Idle
            printing ? println("Ambulance $(SelectedAmb.index) has returned to base at $(SelectedAmb.returnBase) and is going IDLE" ) : nothing
        end
    end
end

#test2=serveEmergence(test,reqArrived)
#reqArrived
Random.seed!(32)
clock  = Clock() 
requestarrival(num_call)  
run!(𝐶, 3000) 

response_times = mean([x.response_time - x.time_request for x in ReqServed])

#evtol que revise en la cola de la región si hay alguien in hard to reach zone y en rango

num_call


dfReqServed = DataFrame(index = [x.index for x in ReqServed], 
                        region = [x.region for x in ReqServed], reqhospit = [x.reqhospit for x in ReqServed],
                        travtimesfrombases = [x.travtimesfrombases for x in ReqServed], travttonearesth = [x.travttonearesth for x in ReqServed],
                        onsitetime = [x.onsitetime for x in ReqServed], traveltfromhtobases = [x.traveltfromhtobases for x in ReqServed],
                        inqueue = [x.inqueue for x in ReqServed], time_request = [x.time_request for x in ReqServed],
                        response_time = [x.response_time for x in ReqServed], 
                        timetoHospital = [x.timetoHospital for x in ReqServed], 
                        finishedService = [x.finishedService for x in ReqServed],
                        servedbyAmbIndex = [x.ServedbyAmbIndex for x in ReqServed])

CSV.write("resultadosF2.csv", dfReqServed) 


