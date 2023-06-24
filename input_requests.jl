module inputreq

include("types.jl")

using .typs
using CSV
using DataFrames

#ReqFileName="test_1_req.csv"

function readReqData2(ReqFileName::String)
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

#ListReq=readReqData2(ReqFileName)

export readReqData2, AmbStatus, Idle, Bussy, Null

end