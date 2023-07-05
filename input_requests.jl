module inputreq

include("types.jl")

using .typs
using CSV
using DataFrames
using DelimitedFiles


ReqFileName="test_1_req.csv"

function readReqData2(ReqFileName)
    dfrequests=DataFrame(CSV.File(ReqFileName))
    travelTimesToSite=[split(x, ",")  for x in dfrequests.travtimesfrombases]
    travtimesfrombases1=[[parse(Float64, s) for s in v] for v in travelTimesToSite]
    travelTimeseReturnBase=[split(x, ",")  for x in dfrequests.traveltfromhtobases]
    travelTimeseReturnBase1=[[parse(Float64, s) for s in v] for v in travelTimeseReturnBase]
    Requests=Array{Request2}(undef,size(dfrequests,1))
    for (i, row) in enumerate(eachrow(dfrequests))
        Requests[i]=Request2(dfrequests[i,:index],
                            dfrequests[i,:region],
                            dfrequests[i,:latitude],
                            dfrequests[i,:longitude],
                            travtimesfrombases1[i],
                            dfrequests[i,:travttonearesth],
                            60.0,
                            travelTimeseReturnBase1[i],
                            dfrequests[i,:priority],
                            dfrequests[i,:reqhospit],
                            dfrequests[i,:inhardtoreachzone],
                            false,0.0,0.0,0.0,0,0.0)
    end
    return Requests
end


function generatearrivals(Requests, freq_calls, freq_horaria_acc, weeks_sim, numReg)
    weekly_calls_reg=reduce(vcat, map(sum, freq_calls))
    sim_reg_calls_reg=weekly_calls_reg.*weeks_sim
    weeks_arr_regs = Vector{Vector{Float64}}(undef, numReg)
    for r in 1:numReg
        weeks_arr_regs_ind = Float64[]
        for j in 1:weeks_sim
            weekly_arr=Float64[]
            for (index, i) in enumerate(freq_calls[r])
                random_numbers = rand(i)
                daily_arr = Float64[]
                for x in random_numbers
                    idx = findfirst(c -> x <= c, freq_horaria_acc)
                    arrival = (((j-1)*24*60*7)+(index-1)*24*60)+round((idx -1)*60 + rand()*60, digits=2)
                    push!(daily_arr, arrival)
                end
                sort!(daily_arr) 
                append!(weekly_arr,daily_arr)
            end
            append!(weeks_arr_regs_ind,weekly_arr)
        end
        weeks_arr_regs[r]=weeks_arr_regs_ind
    end

    arrivals = Vector{Vector{Request2}}(undef, numReg)
    id_arr=0
    for i in 1:numReg
        arrivals[i] = Vector{Request2}()  
        for j in 1:sim_reg_calls_reg[i]
            id_arr+=1
            randomreq=deepcopy(Requests[rand(((1000*i)-999):1000*i)])
            randomreq.index=id_arr
            randomreq.time_request=weeks_arr_regs[i][j]
            push!(arrivals[i],randomreq)
        end
    end

    return arrivals
end

export readReqData2, generatearrivals, AmbStatus, Idle, Bussy, Null

end