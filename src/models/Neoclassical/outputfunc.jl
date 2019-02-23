
function outputfunc(p::NeoClassicalSimple)
    @warn "using deprecated outputfunc NeoClassicalSimple"
    @unpack α = p.params

    vK = p.tStateVectors[1]
    vAlog = p.tStateVectors[2]

    return (vK.^α)*exp.(vAlog')
end
