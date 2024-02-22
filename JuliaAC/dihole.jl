using JuLIP, Plots, LinearAlgebra, NeighbourLists, ASE

include("FIO.jl")
include("AtC.jl")

function plotat(at, i, ii, ibdy, iibdy)
    x, y, z = xyz(at)
    Plots.scatter(x,y,z)
    Plots.scatter!(x[i],y[i],z[i])
    Plots.scatter!(x[ii],y[ii],z[ii])
    Plots.scatter!(x[ibdy],y[ibdy],z[ibdy])
    Plots.scatter!(x[iibdy],y[iibdy],z[iibdy])
end

function plotat(at)
    x, y, z = xyz(at)
    Plots.scatter(x,y,z)
end

function construct_dihole(L, l, lw)
    at = bulk(:W, cubic=true) * L
    # locate the 'center' of atomistic lattice
    Xat = positions(at)
    xcell = diag(cell(at))/2
    xcidx = findall(x->isapprox(x, xcell), Xat)
    if isnothing(xcidx)
        r = [ norm(x - xcell) for x in Xat ]
        _, xcidx = findmin(r)
        Xcidx = xcidx.I[2]
        @assert !isnothing(xcidx)
    end
    xc = Xat[xcidx]
    Xat .-= xc
    set_positions!(at, Xat)

    lp = [0.5 * xcell[1], 0.0, 0.0]
    rp = [-0.5 * xcell[1], 0.0, 0.0]
    r = [ norm(x - lp) for x in Xat ]
    _, il = findmin(r)

    r = [ norm(x - rp) for x in Xat ]
    _, ir = findmin(r)

    nlist = neighbourlist(at, l)
    ildef = nlist.j[findall(nlist.i .== il)]
    irdef = nlist.j[findall(nlist.i .== ir)]

    ildefall = sort([ildef; il])
    irdefall = sort([irdef; ir])
    idef = sort([ildef; il; irdef; ir])

    # inner_nlist = neighbourlist(at, l-4.5, recompute=true)
    inner_nlist = neighbourlist(at, l-lw, recompute=true)
    inner_ildef = inner_nlist.j[findall(inner_nlist.i .== il)]
    inner_irdef = inner_nlist.j[findall(inner_nlist.i .== ir)]

    inner_ildefall = sort([inner_ildef; il])
    inner_irdefall = sort([inner_irdef; ir])
    inner_idef = sort([inner_ildef; il; inner_irdef; ir])

    torus_idef = setdiff(idef, inner_idef)

    at0 = deepcopy(at)

    deleteat!(at, idef)
    return at, at0, idef, ildefall, il, irdefall, ir, inner_idef, inner_ildef, inner_irdef, torus_idef

end

function get_atdef(L, l, lw)
    at, at0, idef, ildef, il, irdef, ir, inner_idef, inner_ildef, inner_irdef, torus_idef = construct_dihole(L, l, lw);
    y0 = copy(at.X)
    Xtorus = at0.X[torus_idef]
    atdef = Atoms(:W, Xtorus)
    return atdef
end

