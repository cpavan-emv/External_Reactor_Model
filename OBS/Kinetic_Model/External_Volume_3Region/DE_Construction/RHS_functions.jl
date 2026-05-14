##############################
#=
The functions in this file are for describing the case of compression with 3 volumes
Volume A is the volume of fixed size and contains the catalyst (reactions allowed)
Volume B is the clearance volume of the engine. Fixed size, no reactions, same P as Volume B
Volume C is the displaced volume. Variable size, no reactions, same P as Volume C
Mass flow from B to A is accounted for on RHS based on finite pressure drop
=#
##############################
function RHS_Common!(du,u,p,t)
    N=Nspec+1
    # p=param=([gasA, gasB, gasC], V_dVfunc, rhocat, [Tw_A, Tw_B, Tw_C], τ)
    # we assume gas state has already been set by u previously
    gasses=p[1]
    V_dV=p[2](t)
    rhocat=p[3]
    Tw=p[4]
    τ=p[5] # rate constant for heat loss (units s)
    XA=[u[3:N]; 1-sum(u[3:N])]
    # assume that catalyst is at wall temperature
    r=[0.0, 0.0, 0.0]
    try
        r.= rhocat*[fluid_props.r1(Tw[1],u[1]/100e3,XA),
            fluid_props.r2(Tw[1],u[1]/100e3,XA), 
            fluid_props.r3(Tw[1],u[1]/100e3,XA)] # mol/s
    catch
    end
    r[isnan.(r)] .= 0.0 # over-ride for case of pure N2 causing NANs
    r[isinf.(r)] .= 0.0 # over-ride for case of pure N2 causing NANs

    # here is all the heat loss terms
    for i=1:3
        k=fluid_props.rho(gasses[i])*fluid_props.cp(gasses[i]) ./ τ[i] # units kW/m3
        du[(i-1)*N+2] = -k*(u[(i-1)*N+2]-Tw[i])
    end

    # determine the upstream direction for compression
    V_dV[2]>0 ? gasup=gasses[2] : gasup=gasses[3]
    rhoC=fluid_props.rho(gasses[3]); Yup=fluid_props.Y(gasup)[spec_ind[1:end-1]]

    # first block is region A -> just reactions
    du[3:N] .= (spec_MW[1:end-1]*1e-3).*[-r[1]+r[2];# production of CO in kg/s
            -r[2]-r[3];# production of CO2
            -2*r[1]-r[2]-3*r[3];# production of H2
            r[1]+r[3];# production of CH3OH
            r[2]+r[3]] # production of H2O,

    # second block is region B
    # only RHS term is the effect of mass flow from C
    du[N+1] += -rhoC*V_dV[2]
    du[N+2] += -fluid_props.enthalpy(gasup)*rhoC*V_dV[2]
    du[N+3:2*N] .+= -rhoC*Yup*V_dV[2]

    # third block is region C
    # region C has exhange with B, no reactions and changing size
    du[2*N+1] = 0 # over-ride this equation for equal pressures
    du[2*N+2] += (fluid_props.enthalpy(gasup)-fluid_props.enthalpy(gasses[3]))*rhoC*V_dV[2]/V_dV[1]
    du[2*N+3:end] .+= (Yup .- fluid_props.Y(gasses[3])[spec_ind[1:end-1]])*rhoC*V_dV[2]/V_dV[1]
end

function f_decoupled!(du, u, p, t)
    Mass=Mass_full(u,p,t)
    RHS = 0.0*du
    RHS_Common!(RHS,u,p,t)
    du .= Mass\RHS
    # if any(isnan.(du))
    #     println(RHS)
    #     println(du)
    #     println(u)
    #     throw(error)
    # end
    return nothing
end

function f_coupled!(du, u, p, t)
    Mass=Mass_full(u,p,t)
    RHS = 0.0*du
    RHS_Common!(RHS,u,p,t)
    N=Nspec+1
    # now add the mass flow across the valves
    # p=param=([gasA, gasB, gasC], V_dVfunc, rhocat, [Tw_A, Tw_B, Tw_C], τ, [K to ext., K from ext., K of intake valve, K of exhaust valve], Vex)
    # direction depends on which pressure is higher
    V_dV=p[2](t)
    Vex=p[7]
    V_dV[2] > 0 ? K=p[6][2] : K=p[6][1]
    mdot_BA = K * (u[N+1] - u[1])
    mdot_BA > 0 ? gasup=p[1][2] : gasup=p[1][1]

    # add to mass conservation equations
    RHS[1] += mdot_BA
    RHS[N+1] -= mdot_BA*Vex/V_dV[1]
    # add to energy conservation equations
    RHS[2] += fluid_props.enthalpy(gasup)*mdot_BA
    RHS[N+2] -= fluid_props.enthalpy(gasup)*mdot_BA*Vex/V_dV[1]
    # add to species conservation equations
    Yup = fluid_props.Y(gasup)[spec_ind[1:end-1]]
    RHS[3:N] .+= mdot_BA*Yup
    RHS[N+3:2*N] .-= mdot_BA*Yup*Vex/V_dV[1]

    du .= Mass\RHS
    return nothing
end

function f_intake_exhaust!(du, u, p, t)
    Mass=Mass_full(u,p,t)
    RHS = 0.0*du
    RHS_Common!(RHS,u,p,t)
    N=Nspec+1
    # now add the mass flow across the valves
    # p=param=([gasA, gasB, gasC], V_dVfunc, rhocat, [Tw_A, Tw_B, Tw_C], τ, [K to ext., K from ext., K of exhaust valve, K of t valve], TPX_exterior)
    # direction depends on which pressure is higher
    V_dV=p[2](t)
    TPX_exterior=p[7]
    # don't need gasA anymore - use it
    gas_ex=p[1][1]
    fluid_props.setTPX(gas_ex,TPX_exterior, spec_ind)
    Pex=TPX_exterior[2]
    V_dV[2] > 0 ? K=p[6][4] : K=p[6][3]
    mdot_in = K * (Pex - u[N+1]) # positive for flow in. Backflow allowed
    mdot_in < 0 ? gasup=p[1][2] : gasup=gas_ex

    # add to mass conservation equations
    RHS[N+1] += mdot_in
    # add to energy conservation equation
    RHS[N+2] += fluid_props.enthalpy(gasup)*mdot_in
    # add to species conservation equations
    Yup = fluid_props.Y(gasup)[spec_ind[1:end-1]]
    RHS[N+3:2*N] .+= mdot_in*Yup

    du .= Mass\RHS
    return nothing
end

function evolve_coupled(u0,param,tsave, condition=nothing)
    prob=ODE.ODEProblem(DE_Model.f_coupled!,u0,(tsave[1],tsave[end]))
    if !isnothing(condition)
        cb=ODE.ContinuousCallback(condition, (integrator)->ODE.terminate!(integrator))
        sol = ODE.solve(prob,ODE.RadauIIA5(autodiff=false),
            p=param, reltol=1e-8, abstol=1e-10, saveat=tsave, callback=cb)
    else
        sol = ODE.solve(prob,ODE.RadauIIA5(autodiff=false),
            p=param, reltol=1e-8, abstol=1e-10, saveat=tsave)
    end
    y=permutedims(stack(sol.u))
    return sol.t, y
end

function evolve_decoupled(u0,param,tsave, condition=nothing)
    prob=ODE.ODEProblem(DE_Model.f_decoupled!,u0,(tsave[1],tsave[end]))
    if !isnothing(condition)
        cb=ODE.ContinuousCallback(condition, (integrator)->ODE.terminate!(integrator))
        sol = ODE.solve(prob,ODE.RadauIIA5(autodiff=false),
            p=param, reltol=1e-8, abstol=1e-10, saveat=tsave, callback=cb)
    else
        sol = ODE.solve(prob,ODE.RadauIIA5(autodiff=false),
            p=param, reltol=1e-8, abstol=1e-10, saveat=tsave)
    end
    y=permutedims(stack(sol.u))
    return sol.t, y
end

function evolve_intake_exhaust(u0,param,tsave)
    prob=ODE.ODEProblem(DE_Model.f_intake_exhaust!,u0,(tsave[1],tsave[end]))
    sol = ODE.solve(prob,ODE.RadauIIA5(autodiff=false),
        p=param, reltol=1e-8, abstol=1e-10, saveat=tsave)
    y=permutedims(stack(sol.u))
    return sol.t, y
end