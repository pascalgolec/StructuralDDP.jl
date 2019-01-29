
# these functions help to initialize the simulation
# they describe the problem the firm faces to choose its initial capital stock
function initializationproblem(p::NeoClassicalSimple, value::Float64, K::Float64)
    @unpack β, C0 = p.params
    r = (1-β)/β
    return value - (1 + r + C0) * K
end
