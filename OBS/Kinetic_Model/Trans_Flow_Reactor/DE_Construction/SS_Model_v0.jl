function flux(rhou, state)
    _,dx=x_dx()
    D=first_order_upwind(dx)
    return D*(rhou.*state)
end

function set_gas_state!(g, u)
    Y=u[4:end]
    h=u[3]
    P=u[1] # Pa
    fluid_props.setHPY(g, (h,P,Y), spec_ind)
    return nothing
end


function source(u, p, g=nothing)
    # for steady-state, u is [H; Y]
    # p[1] is a dictionary where you can pole things like "Pressure"
    du=0*u
    Y=u[2:end]
    h=u[1]
    P=p[:Pressure]
    
    isnothing(g) ? gas=p[:gas] : gas=g
    fluid_props.setHPY(gas, (h,P,Y), spec_ind)
    XA=fluid_props.X(gas)[spec_ind]
    T=fluid_props.T(gas)
    #println((T, P, Y))
    r=[0.0,0.0,0.0]
    try
        r.= p[:rhocat]*[fluid_props.r1(T,P/100e3,XA),
            fluid_props.r2(T,P/100e3,XA), 
            fluid_props.r3(T,P/100e3,XA)] # mol/s
    catch
    end
    r[isnan.(r)] .= 0.0 # over-ride for case of pure N2 causing NANs
    r[isinf.(r)] .= 0.0 # over-ride for case of pure N2 causing NANs
    du[2:end] .+= (spec_MW[1:end-1]*1e-3).*[-r[1]+r[2];# production of CO in kg/s
            -r[2]-r[3];# production of CO2
            -2*r[1]-r[2]-3*r[3];# production of H2
            r[1]+r[3];# production of CH3OH
            r[2]+r[3]] # production of H2O,

    du[1] = -1.0/(25.4e-3*500e-6)*(T-p[:Twall]) # placeholder for temperature loss

    return du

end

function f_SS!(du, u, p)
    # state is saved by cell as default
    # i.e. so that reshape(u, (Nvar,Ncell)) puts the full cell state in each column
    # this is good for the transient side, but is bad for the spatial side
    state_byVar=reshape(u,(NvNc[1],NvNc[2]))'
    state_byCell=reshape(u,(NvNc[1],NvNc[2]))
    p=p[1]
    mdot=p[:mdot_A] * ones(NvNc[2])
    f=stack([flux(mdot, c) for c in eachcol(state_byVar)]) # each column contains the flux for one variable
    # add the inlet flux to the leftmost cell
    _, dx = x_dx()
    f[1,:].-=p[:inlet_state]*p[:mdot_A]/dx[1]
    s=stack([source(c, p, g) for (c, g) in zip(eachcol(state_byCell),p[:gasses])])[:]
    du .=-MapMat'*f[:]+s
    return nothing
end

function j_SS!(J,u,p)
    D=first_order_upwind(x_dx()[2])
    # Jacobian of the differential operator
    Jdiff=SparseArrays.spdiagm(0=>stack(fill(diag(D),NvNc[1]))[:],
        -1=>stack(fill([diag(D,-1);0],NvNc[1]))[1:end-1])
    Jsource=j_source(u,p[1])
    J .= MapMat'*(-Jdiff+Jsource)*MapMat
end

function source_all(u,p)
    return stack([source(c,p) for c in eachcol(reshape(u,(NvNc[1],NvNc[2])))] )[:]
end

function j_source(u0,p, solveP=false)
    #=
    The source term depends only on the local state -> no coupling between cells
    To minimize function calls, I can perturb all cells simultaneously
    =#
    u=copy(u0)
    du=similar(u)
    res=similar(u)
    s0=source_all(u0,p)
    du .= u0*1e-8 .+ 1e-12
    Js=SparseArrays.spzeros(length(u),0)
    for i in 1:NvNc[1]
        u .= copy(u0)
        u[i:NvNc[1]:end] .+= du[i:NvNc[1]:end]
        res.=source_all(u,p)
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
    du .= source(u, p) /p[:mdot_A]
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

function f_SS_P!(du, u, p)
    # same as f_SS but includes momentum equation
    # uses momentum equation to solve for pressure
    # pressure is now the first variable in U
    state_byVar=reshape(u,(NvNc[1],NvNc[2]))'
    state_byCell=reshape(u,(NvNc[1],NvNc[2]))
    p=p[1]
    mdot=p[:mdot_A] * ones(NvNc[2])
    # store the current state varaibles in cantera state
    # this makes accessing them later easier
    [set_gas_state!(g,c) for (g,c) in zip(p[:gasses], eachcol(state_byCell))]
    # each column contains the flux for one variable
    f=[zeros(NvNc[2]) stack([flux(mdot, c) for c in eachcol(state_byVar[:,2:end])])]
    # add the inlet flux to the leftmost cell
    _, dx = x_dx()
    f[1,2:end].-=p[:inlet_state][2:end]*p[:mdot_A]/dx[1]
    s=stack([source_P(g,p) for g in p[:gasses]])[:]
    # add pressure gradient to momentum equation
    # DO NOT DO THIS IN SOURCE TERM -> fucks up the Jacobian calculation
    gradP=central_diff(x_dx()[2][1])*state_byVar[:,1]
    # now stack this correctly
    Ps=[zeros(NvNc[2]);-gradP;state_byVar[:,2].*gradP; zeros((NvNc[1]-3)*NvNc[2])]
    # going to avoid BCs on pressure - I think it should be implicit by the flux BC (we will see...)
    du .=-MapMat'*(f[:]+Ps)+s
    return nothing
end

function source_P(gas, p)
    # this version assumes gas state has already been set
    du=zeros(Nspec+2)
    XA=fluid_props.X(gas)[spec_ind]
    T=fluid_props.T(gas)
    P=fluid_props.P(gas)
    #println((T, P, Y))
    r=[0.0,0.0,0.0]
    #try
        r.= p[:rhocat]*[fluid_props.r1(T,P/100e3,XA),
            fluid_props.r2(T,P/100e3,XA), 
            fluid_props.r3(T,P/100e3,XA)] # mol/s
    #catch
    #end
    #r[isnan.(r)] .= 0.0 # over-ride for case of pure N2 causing NANs
    #r[isinf.(r)] .= 0.0 # over-ride for case of pure N2 causing NANs
    du[4:end] .+= (spec_MW[1:end-1]*1e-3).*[-r[1]+r[2];# production of CO in kg/s
            -r[2]-r[3];# production of CO2
            -2*r[1]-r[2]-3*r[3];# production of H2
            r[1]+r[3];# production of CH3OH
            r[2]+r[3]] # production of H2O,

    du[3] = -1.0/(25.4e-3*500e-6)*(T-p[:Twall]) # placeholder for temperature loss

    return du

end

function j_SS_P!(J,u,p)
    state_byCell=reshape(u,(NvNc[1],NvNc[2]))
    state_byVar=reshape(u,(NvNc[1],NvNc[2]))'
    p=p[1]
    [set_gas_state!(g,c) for (g,c) in zip(p[:gasses], eachcol(state_byCell))]

    D=first_order_upwind(x_dx()[2])
    # Jacobian of the differential operator
    Jdiff=p[:mdot_A]*SparseArrays.spdiagm(0=>stack(fill(diag(D),NvNc[1]))[:],
        -1=>stack(fill([diag(D,-1);0],NvNc[1]))[1:end-1])
    Jsource=j_source_P(u,p)
    #=
    This term needs to look like:
    [0 0 0 0
    -∂x 0 0 0
    u∂x ∂xP 0 0
    0 0 0 0]
    =#

    grad_mat=central_diff(x_dx()[2][1])
    u_grad_mat=diagm(state_byVar[:,2])*grad_mat
    grad_mat_P=grad_mat*state_byVar[:,1]
    Nc=NvNc[2]
    N=NvNc[1]*Nc
    JgradP=SparseArrays.spdiagm(N,N,-Nc=>[Vector(-diag(grad_mat));grad_mat_P],
        -Nc+1=>[0;Vector(-diag(grad_mat,1))],
        -Nc-1=>Vector(-diag(grad_mat,-1)),
        -2*Nc=>diag(u_grad_mat),
        -2*Nc+1=>[0;diag(u_grad_mat,1)],
        -2*Nc-1=>diag(grad_mat,-1))
    

    # right MapMat converts state being acted on to Var-ordered
    # left MapMat converts back to Cell-ordered after operation
    J .= MapMat'*(-Jdiff+JgradP+Jsource)*MapMat 
end

function source_all_P(p)
    return stack([source_P(g,p) for g in p[:gasses]] )[:]
end

function j_source_P(u0,p)
    #=
    The source term depends only on the local state -> no coupling between cells
    To minimize function calls, I can perturb all cells simultaneously
    =#
    u=copy(u0)
    du=similar(u)
    res=similar(u)
    s0=source_all_P(p)
    du .= u0*1e-8 .+ 1e-12
    Js=SparseArrays.spzeros(length(u),0)
    for i in 1:NvNc[1]
        u .= copy(u0)
        u[i:NvNc[1]:end] .+= du[i:NvNc[1]:end]
        state_byCell=reshape(u,(NvNc[1],NvNc[2]))
        [set_gas_state!(g,c) for (g,c) in zip(p[:gasses], eachcol(state_byCell))]
        res.=source_all_P(p)
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
