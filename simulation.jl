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

#GROUND AMBULANCE NETWORK 
AmbFileName="inputAmbulances2.csv" 
dfambulances2=DataFrame(CSV.File(AmbFileName))  #Read the ambulance file
numReg=length(unique(dfambulances2.region)) #Determine the number of regions
Num_estaciones=Vector(combine(groupby(dfambulances2, :region), :station => maximum => :station)[:, "station"])  #Determine the number of stations per region
NumAmbulances=size(dfambulances2)[1] #Determine the number of ambulances in the system
resourcess=inputamb.Inputambs2(dfambulances2, numReg, Num_estaciones, NumAmbulances) #Create an array with the ambulances organized by region, station and type

#Create an array to store the request on queues organized by region and by priority
queuesRegions2=createqueues(numReg)  

# ARRIVALS AND ATRIBUTES OF THE REQUESTS
freq_calls=[[14,  10,  11,  11,  11,  13,  11],  #Weekly frequency of request MTWTFSS
            [9,  6,  7,  7,  7,  8,  7],
            [12,  9,  10,  10,  10,  12,  9],
            [7,  5,  6,  6,  6,  7,  6],
            [18,  13,  15,  15,  15,  17,  14],
            [4,  3,  4,  4,  4,  4,  3],
            [23,  16,  19,  19,  19,  21,  17],
            [53,  38,  44,  44,  44,  50,  41],
            [21,  15,  17,  17,  17,  19,  16],
            [11,  7,  9,  9,  9,  10,  8]]

freq_horaria_acc=[0.027,  0.049,  0.074,  0.096,  0.117,  0.1365, #Cumulative distribution function of probability of receiving a request in each hour of the day
             0.165,  0.196,  0.243,  0.303,  0.365,  0.425,  0.481,  
             0.539,  0.59,  0.638,  0.691,  0.743,  0.797,  0.842,  0.892,  0.932,  0.97,  1]
weeks_sim=8 #Number of weeks in the simulation
ReqFileName="test_1_req.csv"
ListReq=readReqData2(ReqFileName) #Read the potential request data



#Create an array and a function to collect statistics on served requests
ReqServed=Array{Main.inputreq.typs.Request2,1}(undef, 0)
function statisticsReq(ServedReq::Main.inputreq.typs.Request2)
    push!(ReqServed,ServedReq)
end

#Function to select the nearest available eVTOL aircraft considering the range, if the request is in hard-to-reach zone and the region
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
    for k in order_stations 
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

function revisequeues(queuesRegions2,time)
    for reg in queuesRegions2
        #filter!(request -> (time - request.time_request) >= 120, reg[1])
        #filter!(request -> (time - request.time_request) >= 240, reg[2])
        filter!(request -> (time - request.time_request) >= 720, reg[3])
    end
end


function requestarrival(num_call)
    reqArrived=sorted_arrivals[num_call] # selecciona request aleatorio dentro de las llamadas
    global num_call+=1
    event!(𝐶,fun(requestarrival, num_call), at, sorted_arrivals[num_call].time_request) # programar el siguiente request
    revisequeues(queuesRegions2,sorted_arrivals[num_call].time_request)
    order_stations=sortperm(reqArrived.travtimesfrombases)
    if reqArrived.inhardtoreachzone==true && reqArrived.region in regs_eVTOL
        amb_selected=selectNearesteVTOL(ListeVTOL,resourcess,order_stations,reqArrived)
    else
        amb_selected=selectNearestAmbAv(resourcess,order_stations,reqArrived)  #SELECCIONAMOS THE NEAREST AMBULANCE AVAILABLE
    end
    #printing ? println(tau(), ": Request $(reqArrived.index) has arrived and the ambulance selected is $(amb_selected)") : nothing
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
        #printing ? println("Request $(Request.index) has been served by EVTOL: $(SelectedAmb.index)" ) : nothing
        return_to_base=Request.response_time+Request.onsitetime+tt_evtol
        Request.finishedService=return_to_base
        SelectedAmb.returnBase=return_to_base
        statisticsReq(Request)
        event!(𝐶, fun(followingReq,SelectedAmb), at, return_to_base)
    else    
        SelectedAmb.status=Main.inputamb.typs.Bussy
        Request.ServedbyAmbIndex=SelectedAmb.index
        #printing ? println("Request $(Request.index) has been served by Ambulance: $(SelectedAmb.index)" ) : nothing
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
                #printing ? println("EVTOL $(SelectedAmb.index) has returned to base at $(SelectedAmb.returnBase) and is going to: $(followingRequest.index)" ) : nothing
            else
                SelectedAmb.status=Main.auxfunct.typs.Idle     
            end                    
        else
            SelectedAmb.status=Main.auxfunct.typs.Idle
            #printing ? println("EVTOL $(SelectedAmb.index) has returned to base at $(SelectedAmb.returnBase) and is going IDLE" ) : nothing
        end
    else
        if (reduce(+, length.(queuesRegions2[SelectedAmb.region])) > 0)
            followingRequest=priority(SelectedAmb.region,queuesRegions2)
            popfirst!(queuesRegions2[SelectedAmb.region][followingRequest.priority])
            event!(𝐶, fun(serveEmergence,SelectedAmb,followingRequest),at ,SelectedAmb.returnBase)
            #printing ? println("Ambulance $(SelectedAmb.index) has returned to base at $(SelectedAmb.returnBase) and is going to: $(followingRequest.index)" ) : nothing
        else
            SelectedAmb.status=Main.inputamb.typs.Idle
            #printing ? println("Ambulance $(SelectedAmb.index) has returned to base at $(SelectedAmb.returnBase) and is going IDLE" ) : nothing
        end
    end
end


#EVTOL AMBULANCE NETWORK
max_autnonomy=60  #Autonomy of eVTOL aircraft in Km
eVTOL_speed=165 #Speed of eVTOL aircraft in Km/h
LocationsFileName="input_location_samu.csv"
dfLocation=DataFrame(CSV.File(LocationsFileName))  #Create a dataframe with the potential locations for vertiports
n_vertiports=5 #Number of vertiports to be located in the area under study
n_evtols=10  #Number of eVTOL located in the previous vertiports n_evtols>=n_vertiports

################################################# DESDE ACÁ EMPIEZA LA SIMULACION REPLICAS




####################ACA EMPIEZA LA APROXIMACIÓN ALEATORIA
ListeVTOL=randomdecisions(numReg, n_vertiports, n_evtols, Num_estaciones, dfLocation)  #RANDOM DECISION
VectorSol=zeros(Int, size(dfLocation)[1])
for eVTOL in ListeVTOL
    VectorSol[eVTOL.relative_pos]+=1
end 


n_solutions=30 #NUMBER OF DIFFERENT X TESTED
best_replica=Vector{Float64}(undef, n_solutions)
standard_dev=Vector{Float64}(undef, n_solutions)
solutions=Vector{Vector{Main.auxfunct.typs.eVTOL}}(undef, n_solutions)
n_replicas=100
soltime=Vector{Float64}(undef, n_solutions)

for k in 1:n_solutions
    FO=Vector{Float64}(undef, n_replicas)
    clock  = Clock() 
    global num_call=1
    global ListeVTOL=randomdecisions(numReg, n_vertiports, n_evtols, Num_estaciones, dfLocation) #Decisions regarding the location of vertiports and deployment of eVTOL aircraft
    global regs_eVTOL = unique([ev.region for ev in ListeVTOL]) #Determine the regions where the vertiports are located
    ListeVTOL2=deepcopy(ListeVTOL)
    soltime[k]=@elapsed for i in 1:n_replicas
        global num_call=1
        ReqServed=Array{Main.inputreq.typs.Request2,1}(undef, 0)
        queuesRegions2=createqueues(numReg)  
        resourcess=inputamb.Inputambs2(dfambulances2, numReg, Num_estaciones, NumAmbulances)
        ListeVTOL=ListeVTOL2
        #println(ListeVTOL)
        arrivalsF=generatearrivals(ListReq, freq_calls, freq_horaria_acc, weeks_sim, numReg)
        dim1array = reduce(vcat, arrivalsF) 
        global sorted_arrivals = sort(dim1array, by = x -> x.time_request)
        clock  = Clock() 
        requestarrival(num_call)  
        run!(𝐶, 80000) 
        resetClock!(𝐶)
        FO[i] = mean([x.response_time - x.time_request for x in ReqServed])
    end
    best_replica[k]=mean(FO)
    standard_dev[k]=std(FO)
    solutions[k]=ListeVTOL2
end


best_replica
min_value = minimum(best_replica)
mian_value = maximum(best_replica)
using DelimitedFiles
writedlm("resultm4.txt", solutions) 


positions = []
field_region = []
field_station = []

for (i, arr) in enumerate(solutions)
    for x in arr
        push!(positions, i)
        push!(field_region, getproperty(x, :region))
        push!(field_station, getproperty(x, :station))
    end
end

df2= DataFrame(Position = positions, region = field_region, station=field_station)

CSV.write("Solutiones1.csv", df2) 

#HYPERBOX########################

# CONJUNTOS S Y E
Sampling=Vector{Vector{Main.auxfunct.typs.eVTOL}}(undef, 1000)  #SAMPLING VECTOR
Estimator=Vector{Vector{Main.auxfunct.typs.eVTOL}}(undef, 1000) #ESTIMATING VECTOR
a=30  #    SIMULATION EFFORT PER VISITED SOLUTION

#SOLUCION ALEATORIA DADA POR EL USUARIO, ESTO ES PARA X0
ramdom_initial_solution=randomdecisions(numReg, n_vertiports, n_evtols, Num_estaciones, dfLocation) #Decisions regarding the location of vertiports and deployment of eVTOL aircraft
best_sol=Main.auxfunct.typs.eVTOL 
best_sol=ramdom_initial_solution
FO=Vector{Float64}(undef, a)
Sampling[1]=ramdom_initial_solution
for i in 1:a
    global num_call=1
    ReqServed=Array{Main.inputreq.typs.Request2,1}(undef, 0)
    queuesRegions2=createqueues(numReg)  
    resourcess=inputamb.Inputambs2(dfambulances2, numReg, Num_estaciones, NumAmbulances)
    global ListeVTOL=ramdom_initial_solution
    global regs_eVTOL = unique([ev.region for ev in ListeVTOL])
    #println(ListeVTOL)
    arrivalsF=generatearrivals(ListReq, freq_calls, freq_horaria_acc, weeks_sim, numReg)
    dim1array = reduce(vcat, arrivalsF) 
    global sorted_arrivals = sort(dim1array, by = x -> x.time_request)
    clock  = Clock() 
    requestarrival(num_call)  
    run!(𝐶, 80000) 
    resetClock!(𝐶)
    FO[i] = mean([x.response_time - x.time_request for x in ReqServed])
end
G_opt=mean(FO)

nuevoscandidatos=50 
Sampling2=Vector{Vector{Main.auxfunct.typs.eVTOL}}(undef, nuevoscandidatos)
for tm in 1:nuevoscandidatos
    Sampling2[tm]=randomdecisions(numReg, n_vertiports, n_evtols, Num_estaciones, dfLocation) #
end
push!(Sampling2,ramdom_initial_solution)

positionbest=0
Estimator2=unique(Sampling2)  
Y=Vector{Float64}(undef, length(Estimator2))
for (i,tm) in enumerate(Estimator2)
    FO=Vector{Float64}(undef, a)
    for i in 1:a
        global num_call=1
        ReqServed=Array{Main.inputreq.typs.Request2,1}(undef, 0)
        queuesRegions2=createqueues(numReg)  
        resourcess=inputamb.Inputambs2(dfambulances2, numReg, Num_estaciones, NumAmbulances)
        global ListeVTOL=tm
        global regs_eVTOL = unique([ev.region for ev in ListeVTOL])
        arrivalsF=generatearrivals(ListReq, freq_calls, freq_horaria_acc, weeks_sim, numReg)
        dim1array = reduce(vcat, arrivalsF) 
        global sorted_arrivals = sort(dim1array, by = x -> x.time_request)
        clock  = Clock() 
        requestarrival(num_call)  
        run!(𝐶, 80000) 
        resetClock!(𝐶)
        FO[i] = mean([x.response_time - x.time_request for x in ReqServed])
    end
    Y[i]=mean(FO)
    if Y[i]< G_opt
        G_opt=Y[i]
        best_sol=tm
        positionbest=i
    end
end

best_sol
G_opt
positionbest
G_opt_Sol=zeros(Int, size(dfLocation)[1])
for eVTOL in best_sol
    G_opt_Sol[eVTOL.relative_pos]+=1
end 


Estimator2Converted=[[] for _ in 1:length(Estimator2)]
for (i, element) in enumerate(Estimator2)
    VectorSolaux=zeros(Int, size(dfLocation)[1])
    for eVTOL in element
        VectorSolaux[eVTOL.relative_pos]+=1
    end 
    Estimator2Converted[i]=VectorSolaux
end
G_opt_Sol
G_opt_Sol==Estimator2Converted[positionbest]

splice!(Estimator2Converted, positionbest)

#he visitado 51 soluciones
lowerbound= zeros(Int, 42)
upperbound=zeros(Int, 42)

for (d, element) in enumerate(lowerbound)
    elements = [arr[d] for arr in Estimator2Converted]
    upperbound[d]=maximum(elements)
    lowerbound[d]=minimum(elements)
end

#ya tengo el hypercubo, ahora tengo que generar las nuevas soluciones en ese hypercubo, tengo que escoger aleatoriamente M positiones, ver que UP está en esas posiciones
# sumar los UP, si es mayor que el numero de eVTOL, cambiar de posiciones aleatorias, colocar 1 eVTOL en las que tengan cero, el resto asignar aleatoriamente observando el 
# lower bound tambien


using DelimitedFiles
writedlm("vectoroptimo.txt", G_opt_Sol) 
writedlm("lowerbound.txt", lowerbound) 
writedlm("upperbound.txt", upperbound)


sk_mas_1= vector_of_vectors = [Vector{Int}() for _ in 1:nuevoscandidatos]
for nuev in 1:nuevoscandidatos
    while true
        global random_vert2 = randperm(42)[1:n_vertiports]
        if (sum(lowerbound[random_vert2])<=n_evtols) && (sum(upperbound[random_vert2])>=n_evtols)
            break
        end
    end
    lb=lowerbound[random_vert2]
    ub=upperbound[random_vert2]
    random_vert2
    # este vector nos dirá cuantos eVTOL habrán en las ubicaciones de random_vert
    num_evtol_newsol= zeros(Int, n_vertiports)
    # cantidad mínima de eVTOL en esas ubicaciones
    for (j, element) in enumerate(lb)
        if element>1
            num_evtol_newsol[j]=element
        else
            num_evtol_newsol[j]=1
        end
    end
    #crear un vector que me diga cuantos eVTOL caben en cada posición considerando lb y ub
    capacity_evtol_newsol= ub-num_evtol_newsol
    acc_capacity=cumsum(capacity_evtol_newsol)
    positionsnevtol=randperm(sum(capacity_evtol_newsol))[1:n_evtols-sum(num_evtol_newsol)]
    for eVTOL in positionsnevtol
        num_evtol_newsol[searchsortedfirst(acc_capacity, eVTOL)]+=1
    end
    newvector=zeros(Int, size(dfLocation)[1])
    newvector[random_vert2]=num_evtol_newsol
    sk_mas_1[nuev]=newvector
end

Estimator2Converted

#Tenemos que convertir esos Anteriores sk_mas_1 en tipo evtol y correr las simulaciones de nevo

List_evtols2=Vector{Main.auxfunct.typs.eVTOL}()
non_zero_elements = [x for x in Estimator2Converted[1] if x > 0]
non_zero_positions = findall(x -> x > 0, Estimator2Converted[1])
for (k,i) in enumerate(non_zero_positions)
    for j in 1:non_zero_elements[k]
        row=dfLocation[i,:]
        push!(List_evtols2,Main.auxfunct.typs.eVTOL(1001+length(List_evtols2),row.region,row.station,row.idpotentialsmur,Main.auxfunct.typs.Idle,0.0,row.latitude,row.longitude,0.0))
    end
end

Sampling3=Vector{Vector{Main.auxfunct.typs.eVTOL}}(undef, nuevoscandidatos)
for p in 1:nuevoscandidatos
    List_evtols2=Vector{Main.auxfunct.typs.eVTOL}()
    non_zero_elements = [x for x in Estimator2Converted[p] if x > 0]
    non_zero_positions = findall(x -> x > 0, Estimator2Converted[p])
    for (k,i) in enumerate(non_zero_positions)
        for j in 1:non_zero_elements[k]
            row=dfLocation[i,:]
            push!(List_evtols2,Main.auxfunct.typs.eVTOL(1001+length(List_evtols2),row.region,row.station,row.idpotentialsmur,Main.auxfunct.typs.Idle,0.0,row.latitude,row.longitude,0.0))
        end
    end
    Sampling3[p]=List_evtols2
end

Sampling3
Sampling2

push!(Sampling3,best_sol)
Estimator3=unique(Sampling3)

Y=Vector{Float64}(undef, length(Estimator3))
for (i,tm) in enumerate(Estimator3)
    FO=Vector{Float64}(undef, a)
    for i in 1:a
        global num_call=1
        ReqServed=Array{Main.inputreq.typs.Request2,1}(undef, 0)
        queuesRegions2=createqueues(numReg)  
        resourcess=inputamb.Inputambs2(dfambulances2, numReg, Num_estaciones, NumAmbulances)
        global ListeVTOL=tm
        global regs_eVTOL = unique([ev.region for ev in ListeVTOL])
        arrivalsF=generatearrivals(ListReq, freq_calls, freq_horaria_acc, weeks_sim, numReg)
        dim1array = reduce(vcat, arrivalsF) 
        global sorted_arrivals = sort(dim1array, by = x -> x.time_request)
        clock  = Clock() 
        requestarrival(num_call)  
        run!(𝐶, 80000) 
        resetClock!(𝐶)
        FO[i] = mean([x.response_time - x.time_request for x in ReqServed])
    end
    Y[i]=mean(FO)
    if Y[i]< G_opt
        G_opt=Y[i]
        best_sol=tm
        positionbest=i
    end
end

best_sol
G_opt
positionbest



for i in 1:n_replicas
    global num_call=1
    ReqServed=Array{Main.inputreq.typs.Request2,1}(undef, 0)
    queuesRegions2=createqueues(numReg)  
    resourcess=inputamb.Inputambs2(dfambulances2, numReg, Num_estaciones, NumAmbulances)
    ListeVTOL=ListeVTOL2
    println(ListeVTOL)
    arrivalsF=generatearrivals(ListReq, freq_calls, freq_horaria_acc, weeks_sim, numReg)
    dim1array = reduce(vcat, arrivalsF) 
    global sorted_arrivals = sort(dim1array, by = x -> x.time_request)
    clock  = Clock() 
    requestarrival(num_call)  
    run!(𝐶, 80000) 
    resetClock!(𝐶)
    FO[i] = mean([x.response_time - x.time_request for x in ReqServed])
    dfReqServed = DataFrame(index = [x.index for x in ReqServed], 
                       region = [x.region for x in ReqServed], reqhospit = [x.reqhospit for x in ReqServed],
                       travtimesfrombases = [x.travtimesfrombases for x in ReqServed], travttonearesth = [x.travttonearesth for x in ReqServed],
                        onsitetime = [x.onsitetime for x in ReqServed], traveltfromhtobases = [x.traveltfromhtobases for x in ReqServed],
                        inqueue = [x.inqueue for x in ReqServed], time_request = [x.time_request for x in ReqServed],
                        response_time = [x.response_time for x in ReqServed], 
                        timetoHospital = [x.timetoHospital for x in ReqServed], 
                        finishedService = [x.finishedService for x in ReqServed],
                        servedbyAmbIndex = [x.ServedbyAmbIndex for x in ReqServed],
                        responsetime = [(x.response_time -x.time_request) for x in ReqServed],
                        priorityy=[x.priority for x in ReqServed])
    CSV.write("resultadosEVTOL_$i.csv", dfReqServed) 
end


# QUE TENEMOSQUE HACEEER  ENCONTRAR QUE ESTÁ MAL EN LA SIM POR QUE SE VA A LA LUNA
# ESTADÍSTICAS PRELIMINARES

#HYPERBOX










