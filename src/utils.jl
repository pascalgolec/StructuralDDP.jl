inbounds(x::Float64, xmin::Float64, xmax::Float64) = max(xmin, min(x, xmax))

# from DifferentialEquations.jl
macro CSI_str(str)
    return :(string("\x1b[", $(esc(str)), "m"))
end
const TYPE_COLOR = CSI"36"
const NO_COLOR = CSI"0"
