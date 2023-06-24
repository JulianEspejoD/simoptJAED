module auxfunct

include("types.jl")

using .typs

function createqueues(numReg)
    queuesRegions2=Vector{Vector{Vector{Request2}}}(undef, numReg)
    for i in 1:numReg
        queuesRegions2[i] = Vector{Vector{Request2}}(undef, 3)
        for j in 1:3
            queuesRegions2[i][j] = []
        end
    end
    return queuesRegions2
end 

export createqueues
end