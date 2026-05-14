function calc_mflux(RHS,p,t)
    dρ_rxn=drho_all(p[:gasses])*RHS # should be a vector of length Ncell
    # ^^ Only included effects of changing temperature and composition
    # now we construct a linear system
    #=
    dρ/dt=-d/dx(rho*u)
    we have Ncell equations dρ=F_i+1/2 - F_i-1/2
    plus a global mass conservation equation (dM/dt = -F_N+1/2 + F_N-1/)
    giving Ncell+1 equations for Ncell+1 face fluxes
    BC is mass flow rate in

    dρ calculated above DOES NOT include pressure contribution
    We need to add, to dρ for each cell, ρ/T * dP
    Either prescribe dP, (then global mass conservation not really needed)
    Or prescribe mdot out (then dM/dt = sum(dρ/dt * Δx)=sum((dρ/dt_rxn + dρ/dt_P)*Δx ))
    =#
    x,dx=x_dx()
    # this matrix is [BC equation; F_i+1/2 - F_i-1/2]
    # and solves for all boundary fluxes
    Mat=SparseArrays.diagm(0 => ones(NvNc[2]+1),
        -1 => -ones(NvNc[2]))
    s = [p[:mdot_A]; -dρ_rxn .* dx]
    if :dP_dt in keys(p)
        s[2:end] .+= [fluid_props.rho(g)/fluid_props.P(g) for g in p[:gasses]]*p[:dP_dt](t) .* dx
    elseif :mdot_A_out in keys(p)
        error("Pressure state variable not implemented")
    else
        error("Missing BC on flux")
    end
    return Mat \ s

end

function f_trans_simple!(du, u, p,t)
    # identical to f_SS but divides through by density
    state_byVar=reshape(u,(NvNc[1],NvNc[2]))'
    state_byCell=reshape(u,(NvNc[1],NvNc[2]))
    p=p[1]
    [set_gas_state!(g,c, p[:Pressure]) for (g,c) in zip(p[:gasses], eachcol(state_byCell))]
    ρ=[fluid_props.rho(g) for g in p[:gasses]]
    mdot=p[:mdot_A] * ones(NvNc[2])
    u=mdot ./ ρ
    # calculate speed:
    f=stack([u .* flux(1.0, c) for c in eachcol(state_byVar)]) # each column contains the flux for one variable
    # add the inlet flux to the leftmost cell
    _, dx = x_dx()
    f[1,:].-=p[:inlet_state]*p[:mdot_A]/(dx[1]*ρ[1])

    # use "Pressure" for a constant pressure and "P" for a functino
    :P in keys(p) ? P=p[:P](t) : P=p[:Pressure]

    s=stack([source(p, g, P) ./ fluid_props.rho(g)  for g in p[:gasses]])[:]
    # mass matrix as constructed applies to cell-ordered
    # multiply by map mat on right to convert arg to var-ordered
    # multiply by map mat transpose on left to convert result back to cell-ordered
    # same thing I do for the jacobian
    #Mass_inv=MapMat'*SparseArrays.spdiagm(repeat(1.0 ./ ρ, NvNc[1]))*MapMat

    #du .=Mass_inv*(-MapMat'*f[:]+s)
    du .=-MapMat'*f[:]+s
    return nothing
end


function f_trans!(du, u, p,t)
    # state is saved by cell as default
    # i.e. so that reshape(u, (Nvar,Ncell)) puts the full cell state in each column
    # this is good for the transient side, but is bad for the spatial side
    state_byVar=reshape(u,(NvNc[1],NvNc[2]))'
    state_byCell=reshape(u,(NvNc[1],NvNc[2]))
    p=p[1]
    [set_gas_state!(g,c, p[:Pressure]) for (g,c) in zip(p[:gasses], eachcol(state_byCell))]
    ρ=[fluid_props.rho(g) for g in p[:gasses]]
    # first step -> calculate source term
    :P in keys(p) ? P=p[:P](t) : P=p[:Pressure]
    s=stack([source(p, g, P)/fluid_props.rho(g) for g in p[:gasses]])[:]
    _, dx = x_dx()
    # second step -> use this to get predictor mass flux
    # these mass fluxes will be on cell faces
    for i in 1:2
        if i>1
            # subsequent passes, update as needed
            mf_face=calc_mflux(du, p,t)
            mf_center=0.5*(mf_face[1:end-1] + mf_face[2:end])
            u=mf_center ./ ρ
        else
            # first pass use the SS value
            u = p[:mdot_A] ./ ρ
        end
        # each column contains the flux for one variable
        f=stack([u .* flux(1.0, c) for c in eachcol(state_byVar)]) 
        # add the inlet flux to the leftmost cell
        f[1,:].-=p[:inlet_state]*p[:mdot_A]/(dx[1]*ρ[1])
        du .=-MapMat'*f[:]+s
    end
    return nothing
end

