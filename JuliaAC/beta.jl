
# maybe compute once before the minimization...
function eval_blending(r::Array{Float64,1}, R1::Float64, R2::Float64, flag::Symbol)
    β = zeros(length(r))
    if flag == :none
        fill!(β, 1.0)
    elseif flag == :shock
        β[findall(r .> R2)] .= 1.0
    elseif flag == :affine
        β[findall(r .> R2)] .= 1.0
        ib = findall(R1 .< r .< R2)
        β[ib] = (r[ib].-R1)/(R2-R1)
    elseif flag == :cos
        β[findall(r .> R2)] .= 1.0
        ib = findall(R1 .< r .< R2)
        β[ib] = cos.(2*(r[ib].-R1)/(π*(R2-R1)))
    else
        @info(`Not implemented yet.`)   
    end
    return β
end

