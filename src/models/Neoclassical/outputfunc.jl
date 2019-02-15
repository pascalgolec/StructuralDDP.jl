
function outputfunc(p::NeoClassicalSimple)
    @unpack α = p.params

    vK = p.tStateVectors[1]
    vAlog = p.tStateVectors[2]

    return (vK.^α)*exp.(vAlog')
end

function outputfunc(p::NeoClassicalSimple)

    rewardfunc(p::DDM, vStateVars::Vector{Float64}, choice)

    @unpack α = p.params

    vK = p.tStateVectors[1]
    vAlog = p.tStateVectors[2]

    return (vK.^α)*exp.(vAlog')
end
