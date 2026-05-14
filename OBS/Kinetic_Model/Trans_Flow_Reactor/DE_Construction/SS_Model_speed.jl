function flux(rhou, state)
    _,dx=x_dx()
    D=first_order_upwind(dx)
    return D*(rhou.*state)
end

function set_gas_state!(g, u, P)
    Y=u[3:end]
    h=u[2]
    fluid_props.setHPY(g, (h,P,Y), spec_ind)
    return nothing
end

# this is point-wise
function drho(g)
    D=LinearAlgebra.diagm(ones(Nspec))
    FB=LinearAlgebra.diagm(ones(Nspec))
    D_hY!(D,g)
    FB_hY!(FB, g)
    return D * LinearAlgebra.inv(FB)
end

# this is point-wise in-place
function drho!(dρ, g)
    D=LinearAlgebra.diagm(ones(Nspec))
    FB=LinearAlgebra.diagm(ones(Nspec))
    D_hY!(D,g)
    FB_hY!(FB, g)
    dρ .= D * LinearAlgebra.inv(FB)
end

function du_const_mdot_point(g, mdot_A)
    return [0.0;drho(g) * ones(Float64, Nspec)]
end

function Jac_drho_u(u, p)
    # assume that the state has already been set
    # each column contains the flux for one variable. returning using transpose makes it so columns are cells
    speed=reshape(u, (NvNc...))'[:] # this is a vector that is variable-ordered
    speed[NvNc[2]+1:end] .= 0.0 # for this function, only want speed
    # Construct an array that is dρ/dϕ = [[dρ_1/dϕ], [dρ_2/dϕ], ...]
    R=collect(1:NvNc[1]*NvNc[2])
    C=stack([i*ones(Int,NvNc[1]) for i in 1:NvNc[2]])[:]
    V=stack([du_const_mdot_point(g, p[:mdot_A]) for g in p[:gasses]])[:]
    tmp=SparseArrays.sparse(R, C, V)
    #mat=hcat(fill(tmp, NvNc[1])...)*MapMat # repeat this for each variable & convert shape to cell-ordered
    mat=hcat(fill(tmp, NvNc[1])...) # repeat this for each variable    
    ρ=[fluid_props.rho(g) for g in p[:gasses]]
    return mat' .* speed + (MapMat'*LinearAlgebra.diagm([ρ; zeros(NvNc[2]*(NvNc[1]-1))]))' # add the diagonal part
end

function source(p, gas)
    # this version assumes gas state has already been set
    du=zeros(Nspec+1)
    XA=fluid_props.X(gas)[spec_ind]
    T=fluid_props.T(gas)
    P=p[:Pressure]
    r=[0.0,0.0,0.0]
    #try
        r.= p[:rhocat]*[fluid_props.r1(T,P/100e3,XA),
            fluid_props.r2(T,P/100e3,XA), 
            fluid_props.r3(T,P/100e3,XA)] # mol/s
    #catch
    #end
    #r[isnan.(r)] .= 0.0 # over-ride for case of pure N2 causing NANs
    #r[isinf.(r)] .= 0.0 # over-ride for case of pure N2 causing NANs
    du[3:end] .+= (spec_MW[1:end-1]*1e-3).*[-r[1]+r[2];# production of CO in kg/s
            -r[2]-r[3];# production of CO2
            -2*r[1]-r[2]-3*r[3];# production of H2
            r[1]+r[3];# production of CH3OH
            r[2]+r[3]] # production of H2O,

    du[2] = -1.0/(25.4e-3*500e-6)*(T-p[:Twall]) # placeholder for temperature loss
    #return du/fluid_props.rho(gas)
    return du
end

function f_SS!(du, u, p)
    # state is saved by cell as default
    # i.e. so that reshape(u, (Nvar,Ncell)) puts the full cell state in each column
    # this is good for the transient side, but is bad for the spatial side
    state_byVar=reshape(u,(NvNc[1],NvNc[2]))'
    state_byCell=reshape(u,(NvNc[1],NvNc[2]))
    p=p[1]
    [set_gas_state!(g,c, p[:Pressure]) for (g,c) in zip(p[:gasses], eachcol(state_byCell))]
    ρ=[fluid_props.rho(g) for g in p[:gasses]]
    mdot=p[:mdot_A] * ones(NvNc[2])
    # calculate speed:
    f=stack([flux(mdot, c) for c in eachcol(state_byVar)]) # each column contains the flux for one variable
    # over-ride first equation -> that is mdot=rho*u
    f[:,1]=mdot - state_byVar[:,1].*ρ
    # add the inlet flux to the leftmost cell
    _, dx = x_dx()
    f[1,2:end].-=p[:inlet_state]*p[:mdot_A]/(dx[1])#*ρ[1])
    s=stack([source(p, g) for g in p[:gasses]])[:]
    du .=-MapMat'*f[:]+s
    return nothing
end

function j_SS!(J,u,p)
    p=p[1]
    #speed=p[:mdot_A] ./ [fluid_props.rho(g) for g in p[:gasses]]
    [set_gas_state!(g,c, p[:Pressure]) for (g,c) in zip(p[:gasses], eachcol(reshape(u,(NvNc...))))]
    D=p[:mdot_A] .*first_order_upwind(x_dx()[2])
    # Jacobian of the differential operator
    Jdiff=SparseArrays.spdiagm(0=>stack(fill(diag(D),NvNc[1]))[:],
        -1=>stack(fill([diag(D,-1);0],NvNc[1]))[1:end-1])
    Jsource=j_source(u,p)
    Ju=Jac_drho_u(u,p)
    J .= MapMat'*(-Jdiff+Jsource)*MapMat
    J[1:NvNc[2],:] .=Ju[1:NvNc[2],:]

end

function source_all(p)
    return stack([source(p,g) for g in p[:gasses]])[:]
end

function j_source(u0,p)
    #=
    The source term depends only on the local state -> no coupling between cells
    To minimize function calls, I can perturb all cells simultaneously
    =#
    u=copy(u0)
    du=similar(u)
    res=similar(u)
    s0=source_all(p)
    du .= u0*1e-8 .+ 1e-12
    Js=SparseArrays.spzeros(length(u),0)
    for i in 1:NvNc[1]
        u .= copy(u0)
        u[i:NvNc[1]:end] .+= du[i:NvNc[1]:end]
        [set_gas_state!(g,c, p[:Pressure]) for (g,c) in zip(p[:gasses], eachcol(reshape(u,(NvNc...))))]
        res.=source_all(p)
        res .-= s0
        # res is ordered by cell, so:
        for j in 1:NvNc[2]
            res[(j-1)*NvNc[1]+1:j*NvNc[1]] .*= 1 ./ du[i:NvNc[1]:end][j]
        end
        # for state vector ordered by variable, we get nice block diagonals
        res_byVar=reshape(res,(NvNc[1],NvNc[2]))' # each column is one variable
        Js=[Js SparseArrays.spdiagm(length(res),NvNc[2],
            ([ind=>r for (ind, r) in zip(-(0:NvNc[2]:length(res)),eachcol(res_byVar))])...)] # diagonal for the ith variable
    end
    return Js
end

function f_SS_IVP!(du, u, p, t)
    p=p[1]
    gas=p[:gasses][1] # only using first gas for IVP
    set_gas_state!(gas, [0;u], p[:Pressure])
    du .= source(p,gas)[2:end] /p[:mdot_A]
    return nothing 
end

function jacobian_brute_force(J,u,p, fun!)
    du=similar(u)
    res=similar(u)
    fun!(res, u, p)
    res0=copy(res)
    du .= u*1e-8 .+ 1e-12
    for i in 1:length(u)
        u_pert = copy(u)
        u_pert[i] += du[i]
        fun!(res, u_pert, p)
        res .-= res0
        res[abs.(res).<eps()] .= 0.0
        res .*= 1 ./ du[i]
        J[:,i] .= res
    end
    return nothing
end
