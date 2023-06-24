
using DiscreteEvents, Random, Distributions
using CSV
using DataFrames
using DataStructures
using StatsBase


#include("definitions.jl")
#include("types.jl")
include("input_ambulances.jl")
include("input_requests.jl")
include("auxfunctions.jl")
using .inputamb
using .inputreq
using .auxfunct

AmbFileName="inputAmbulances2.csv"
resourcessLis=Inputambs(AmbFileName)
NumReg=numRegiones(AmbFileName)
queuesRegions2=createqueues(NumReg)


ReqFileName="test_1_req.csv"
ListReq=readReqData2(ReqFileName)

printing = true
num_call=1

#ACA VAMOS A CREAR UN ARREGLO PARA COLECCIONAR ESTADÍSTICAS DEL SISTEMA

ReqServed=Array{Request2,1}(undef, 0)
function statisticsReq(ServedReq::Request2)
    push!(ReqServed,ServedReq)
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


indexreqArrived=rand(1:length(ListReq))
reqArrived=ListReq[indexreqArrived]
order_stations=sortperm(reqArrived.travtimesfrombases)

test=selectNearestAmbAv(resourcess,order_stations,reqArrived)

function requestarrival(μ)
    indexreqArrived=rand(1:length(ListReq))
    reqArrived=ListReq[indexreqArrived] # selecciona request aleatorio dentro de las llamadas
    deleteat!(ListReq, indexreqArrived)
    ta = rand(Erlang())*μ
    #push!(queuesRegions[(reqArrived.region)],reqArrived)# Coloca en la cola el request
    global num_call+=1
    reqArrived.time_request=tau()
    event!(𝐶,fun(requestarrival, μ), after, ta) # programar el siguiente request
    order_stations=sortperm(reqArrived.travtimesfrombases)
    amb_selected=selectNearestAmbAv(resourcess,order_stations,reqArrived)  #SELECCIONAMOS THE NEAREST AMBULANCE AVAILABLE
    printing ? println(tau(), ": Request $(reqArrived.index) has arrived and the ambulance selected is $(amb_selected)") : nothing
    if amb_selected===nothing
        reqArrived.inqueue=true
        push!(queuesRegions2[reqArrived.region][reqArrived.priority],reqArrived)  #Coloca en la cola el request
    else
        serveEmergence(amb_selected,reqArrived)
    end
end


function serveEmergence(SelectedAmb,Request)# A ESTA FUNCIÓN LE MANDAMOS LA AMBULANCIA SELECCIONADA Y EL TRABAJO
    SelectedAmb.status=Main.inputamb.typs.Bussy
    Request.ServedbyAmbIndex=SelectedAmb.index
    printing ? println("Request $(Request.index) has been served by Ambulance: $(SelectedAmb.index)" ) : nothing
    if (Request.inqueue==false)
        Request.response_time=Request.travtimesfrombases[SelectedAmb.station]+Request.time_request #calculo de tiempo de respuesta aca llegamos al sitio
    else
        Request.response_time=Request.travtimesfrombases[SelectedAmb.station]+SelectedAmb.returnBase
    end
    #statisticsReq(Request)
    #printing ? println("Request $(Request.index) has been served at: $(Request.response_time) by Ambulance: $(SelectedAmb.index)" ) : nothing
    #aca revisamos si tenemos que llevarlo o no al hospital, si es tipo 1 puede si es 2 no
    if  (SelectedAmb.class==2) 
        return_to_base=Request.response_time+Request.onsitetime+Request.traveltfromhtobases[SelectedAmb.station] # TIEMPO DE RETORNO A BASE (response time mas + onsite time + devolucion a base)
        SelectedAmb.returnBase=return_to_base
        Request.finishedService=return_to_base
        #statisticsReq(Request)
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
            #statisticsReq(Request)
            event!(𝐶, fun(followingReq,SelectedAmb), at, return_to_base) # ACA PROGRAMAMOS EL EVENTO CUANDO VUELVA A BASE
            #printing ? println(return_to_base, ": Ambulance $(SelectedAmb.index) has returned to base") : nothing
            return nothing
        end
    end
end


function priority(region,queuesRegions2)
    if !isempty(queuesRegions2[region][1])
        return queuesRegions2[region][1][1]
    elseif  !isempty(queuesRegions2[region][2]) 
        return queuesRegions2[region][2][1]
    else 
        return queuesRegions2[region][3][1]     
    end
end

function followingReq(SelectedAmb)# #aca establece política de prioridad!!! EN ESTA FUNCION CAMBIAMOS EL STATUS DE LA AMBULANCIA 
    if (reduce(+, length.(queuesRegions2[SelectedAmb.region])) > 0)
        followingRequest=priority(SelectedAmb.region,queuesRegions2)
        popfirst!(queuesRegions2[SelectedAmb.region][followingRequest.priority])
        event!(𝐶, fun(serveEmergence,SelectedAmb,followingRequest),at ,SelectedAmb.returnBase)
        printing ? println("Ambulance $(SelectedAmb.index) has returned to base at $(SelectedAmb.returnBase) and is going to: $(followingRequest.index)" ) : nothing
    else
        SelectedAmb.status=Idle
        printing ? println("Ambulance $(SelectedAmb.index) has returned to base at $(SelectedAmb.returnBase) and is going IDLE" ) : nothing
    end
end

test2=serveEmergence(test,reqArrived)
reqArrived
Random.seed!(220)
requestarrival(14)  
run!(𝐶, 3000) 