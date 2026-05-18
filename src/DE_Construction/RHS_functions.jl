##############################
#=
Three-volume compression model:
  Volume A -- fixed external volume, contains catalyst (reactions allowed)
  Volume B -- clearance volume, fixed size, no reactions, same P as C
  Volume C -- displaced variable volume, no reactions, same P as B
  Mass flow B <-> A accounted for on RHS via finite pressure drop
=#
##############################

# Stoichiometry for the five reaction species (CO, CO2, H2, CH3OH, H2O) across
# the three reactions (r1=CO-hyd, r2=RWGS, r3=CO2-hyd). Rows = rxn species,
# cols = reactions. Stored as a tuple-of-tuples so it lives in static memory.
#   r1: CO + 2H2 -> CH3OH
#   r2: CO2 + H2 -> CO + H2O
#   r3: CO2 + 3H2 -> CH3OH + H2O

const _N2_=6950 # q in sLm, P in bar, T in K
const _rho_ref = 1.18 # kg/m3 -> reference density

function sigmoid(x, x0)
    return 0
end
# from swagelok sizing buletin
# https://www.swagelok.com/downloads/webcatalogs/EN/MS-06-84.pdf
function _mdot_cv(Cv, P_up, P_dn, T_up, M_up)
    x =(P_up - P_dn) / P_up
    deltaP_bar=(P_up-P_dn)/100e3
    rhog=P_up/(T_up*ct.Ru/M_up)
    Gg=rhog/_rho_ref
    #if x >= 0.5
        q1= _N2_ * Cv * P_up/100e3 * (1-2*deltaP_bar/(3*P_up/100e3))*
            sqrt(deltaP_bar/(P_up/100e3*Gg*T_up))
    #else
        q2=0.471*_N2_ *Cv*P_up/100e3*sqrt(1/(Gg*T_up))
    #end

    q=q1*0.5*(tanh.(-(x - 0.5)*50)+1) + q2*0.5*(tanh((x-0.5)*50)+1) 

    return q/(60*1000)*101.325e3/(ct.Ru*298.15)*M_up 
end

const _rxn_stoich = ((-1, 1, 0),    # CO
                     (0, -1, -1),   # CO2
                     (-2, -1, -3),  # H2
                     (1, 0, 1),     # CH3OH
                     (0, 1, 1))     # H2O

function RHS_Common!(RHS, u, params::ReactorParams, t)
    gases = params.gases
    V_dV   = params.Vfunc(t)
    Nspec  = gases[1].gas.Nspec
    N      = Nspec + 1
    rhocat = params.rhocat
    Tw     = params.T_walls
    tau    = params.tau
    r      = params.scratch.r

    # Kinetics in region A -- catalyst at wall temperature, composition from current gas state
    r .= rhocat .* [ct.MeOH_kinetics.r1(Tw[1], u[1]/100e3, gases[1].X),
                    ct.MeOH_kinetics.r2(Tw[1], u[1]/100e3, gases[1].X),
                    ct.MeOH_kinetics.r3(Tw[1], u[1]/100e3, gases[1].X)]
    r[isnan.(r)] .= 0.0  # pure-N2 guard
    r[isinf.(r)] .= 0.0

    # Heat loss to walls for all three regions
    for i in 1:2
        k = gases[i].rho * gases[i].cp / tau[i]
        RHS[(i-1)*N+2] = -k * (u[(i-1)*N+2] - Tw[i])
    end

    # Region A: species production. rxn_spec_ind[i] gives the 1-based mechanism position of
    # each reaction species, which equals its position in the state-vector species block.
    # Bath gas (mechanism species Nspec) is not in the state vector -- skip it.
    for i in 1:5
        j = ct.MeOH_kinetics.rxn_spec_ind[i]
        j == Nspec && continue
        nu_i = _rxn_stoich[i][1]*r[1] + _rxn_stoich[i][2]*r[2] + _rxn_stoich[i][3]*r[3]
        RHS[2+j] += gases[1].MW_spec[j] * 1e-3 * nu_i*params.beta
    end

    # Region B: changing volume, no reactions
    RHS[N+1]   -= gases[2].rho * V_dV[2]
    RHS[N+2]   -= gases[2].enthalpy * gases[2].rho * V_dV[2]
    RHS[N+3:end] .-= gases[2].Y[1:end-1]*gases[2].rho * V_dV[2]
end

function f_decoupled!(du, u, p, t)
    params, _ = p
    Mass = Mass_full(u, params, t)
    RHS  = params.scratch.RHS
    fill!(RHS, 0.0)
    RHS_Common!(RHS, u, params, t)
    du .= Mass \ RHS
    return nothing
end

function f_coupled!(du, u, p, t)
    params, _ = p
    Mass = Mass_full(u, params, t)
    RHS  = params.scratch.RHS
    gases=params.gases
    fill!(RHS, 0.0)
    RHS_Common!(RHS, u, params, t)
    N = params.gases[1].gas.Nspec + 1

    # Mass flow between B and A via Cv valve model
    # Cv_vals[1]: reactor-inlet valve (compression stroke, dV/dt < 0)
    # Cv_vals[2]: reactor-outlet valve (expansion stroke, dV/dt > 0)
    V_dV = params.Vfunc(t)
    Cv   = V_dV[2] > 0 ? params.Cv_vals[2] : params.Cv_vals[1]
    P_A, P_B = u[1], u[N+1]
    if P_B >= P_A
        mdot_BA =  _mdot_cv(Cv, P_B, P_A, u[N+2], gases[2].MW_mix)
    else
        mdot_BA = -_mdot_cv(Cv, P_A, P_B, u[2],   gases[1].MW_mix)
    end
    gasup = mdot_BA > 0 ? gases[2] : gases[1]
    mdot_BA /= params.Vd

    RHS[1]      += mdot_BA
    RHS[N+1]    -= mdot_BA
    RHS[2]      += gasup.enthalpy * mdot_BA
    RHS[N+2]    -= gasup.enthalpy * mdot_BA
    Yup          = gasup.Y[1:end-1]
    RHS[3:N]    .+= mdot_BA * Yup
    RHS[N+3:2*N] .-= mdot_BA * Yup

    du .= Mass \ RHS
    return nothing
end

function f_intake_exhaust!(du, u, p, t)
    params, mode = p
    Mass = Mass_full(u, params, t)
    RHS  = params.scratch.RHS
    gases=params.gases
    fill!(RHS, 0.0)
    RHS_Common!(RHS, u, params, t)
    N = params.gases[1].gas.Nspec + 1

    V_dV         = params.Vfunc(t)
    TPX_exterior = mode.TPX_exterior

    # Use gas A slot as scratch for exterior state.
    # TPX_exterior must carry the full Nspec-element X vector.
    # Cv_vals[3]: exhaust valve (compression stroke, dV/dt < 0)
    # Cv_vals[4]: intake valve  (expansion stroke, dV/dt > 0)
    gas_ex = params.gases[1]
    ct.setTPX(gas_ex, TPX_exterior)
    Pex = TPX_exterior[2]
    Cv  = V_dV[2] > 0 ? params.Cv_vals[4] : params.Cv_vals[3]
    P_B = u[N+1]
    if Pex >= P_B
        mdot_in =  _mdot_cv(Cv, Pex, P_B, TPX_exterior[1], gas_ex.MW_mix)
    else
        mdot_in = -_mdot_cv(Cv, P_B, Pex, u[N+2], gases[2].MW_mix)
    end
    gasup = mdot_in < 0 ? gases[2] : gas_ex
    mdot_in /= params.Vd

    RHS[N+1]     += mdot_in
    RHS[N+2]     += gasup.enthalpy * mdot_in
    Yup           = gasup.Y[1:end-1]
    RHS[N+3:2*N] .+= mdot_in * Yup

    du .= Mass \ RHS
    return nothing
end

# Unified evolve with multiple dispatch on ReactorMode

function _run_ode(prob, p, tsave, condition, alg)
    if !isnothing(condition)
        cb  = ODE.ContinuousCallback(condition, integrator -> terminate!(integrator))
        sol = ODE.solve(prob, alg, p=p, reltol=1e-10, abstol=1e-12, saveat=tsave, callback=cb)
    else
        sol = ODE.solve(prob, alg, p=p, reltol=1e-10, abstol=1e-12, saveat=tsave)
    end
    return sol.t, permutedims(stack(sol.u))
end

function evolve(::Decoupled, u0, params::ReactorParams, tsave, condition=nothing;
                alg=_default_alg)
    p    = (params, Decoupled())
    prob = ODE.ODEProblem(f_decoupled!, u0, (tsave[1], tsave[end]))
    return _run_ode(prob, p, tsave, condition, alg)
end

function evolve(mode::Coupled, u0, params::ReactorParams, tsave, condition=nothing;
                alg=_default_alg)
    p    = (params, mode)
    prob = ODE.ODEProblem(f_coupled!, u0, (tsave[1], tsave[end]))
    return _run_ode(prob, p, tsave, condition, alg)
end

function evolve(mode::IntakeExhaust, u0, params::ReactorParams, tsave;
                alg=_default_alg)
    p    = (params, mode)
    prob = ODE.ODEProblem(f_intake_exhaust!, u0, (tsave[1], tsave[end]))
    return _run_ode(prob, p, tsave, nothing, alg)
end
