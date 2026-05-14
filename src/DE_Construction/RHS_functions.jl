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
const _rxn_stoich = ((-1, 1, 0),    # CO
                     (0, -1, -1),   # CO2
                     (-2, -1, -3),  # H2
                     (1, 0, 1),     # CH3OH
                     (0, 1, 1))     # H2O

function RHS_Common!(RHS, u, params::ReactorParams, t)
    gasses = params.gasses
    V_dV   = params.Vfunc(t)
    Nspec  = gasses[1].gas.Nspec
    N      = Nspec + 1
    rhocat = params.rhocat
    Tw     = params.T_walls
    tau    = params.tau
    r      = params.scratch.r

    # Kinetics in region A -- catalyst at wall temperature, composition from current gas state
    r .= rhocat .* [ct.MeOH_kinetics.r1(Tw[1], u[1]/100e3, gasses[1].X),
                    ct.MeOH_kinetics.r2(Tw[1], u[1]/100e3, gasses[1].X),
                    ct.MeOH_kinetics.r3(Tw[1], u[1]/100e3, gasses[1].X)]
    r[isnan.(r)] .= 0.0  # pure-N2 guard
    r[isinf.(r)] .= 0.0

    # Heat loss to walls for all three regions
    for i in 1:3
        k = gasses[i].rho * gasses[i].cp / tau[i]
        RHS[(i-1)*N+2] = -k * (u[(i-1)*N+2] - Tw[i])
    end

    # Upstream direction for compression/expansion flow between B and C
    gasup = V_dV[2] > 0 ? gasses[2] : gasses[3]
    rhoC  = gasses[3].rho
    Yup   = gasup.Y[1:end-1]

    # Region A: species production. rxn_spec_ind[i] gives the 1-based mechanism position of
    # each reaction species, which equals its position in the state-vector species block.
    # Bath gas (mechanism species Nspec) is not in the state vector -- skip it.
    for i in 1:5
        j = ct.MeOH_kinetics.rxn_spec_ind[i]
        j == Nspec && continue
        nu_i = _rxn_stoich[i][1]*r[1] + _rxn_stoich[i][2]*r[2] + _rxn_stoich[i][3]*r[3]
        RHS[2+j] += gasses[1].MW_spec[j] * 1e-3 * nu_i
    end

    # Region B: mass exchange with C
    RHS[N+1]     += -rhoC * V_dV[2]
    RHS[N+2]     += -gasup.enthalpy * rhoC * V_dV[2]
    RHS[N+3:2*N] .+= -rhoC * Yup * V_dV[2]

    # Region C: exchange with B, changing volume (no reactions)
    RHS[2*N+1]    = 0.0  # pressure equation replaced by equal-pressure constraint
    RHS[2*N+2]   += (gasup.enthalpy - gasses[3].enthalpy) * rhoC * V_dV[2] / V_dV[1]
    RHS[2*N+3:end] .+= (Yup .- gasses[3].Y[1:end-1]) * rhoC * V_dV[2] / V_dV[1]
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
    params, mode = p
    Mass = Mass_full(u, params, t)
    RHS  = params.scratch.RHS
    fill!(RHS, 0.0)
    RHS_Common!(RHS, u, params, t)
    N = params.gasses[1].gas.Nspec + 1

    # Mass flow between B and A via finite pressure drop
    V_dV = params.Vfunc(t)
    Vex  = mode.Vex
    K    = V_dV[2] > 0 ? params.K_vals[2] : params.K_vals[1]
    mdot_BA = K * (u[N+1] - u[1])
    gasup   = mdot_BA > 0 ? params.gasses[2] : params.gasses[1]

    RHS[1]      += mdot_BA
    RHS[N+1]    -= mdot_BA * Vex / V_dV[1]
    RHS[2]      += gasup.enthalpy * mdot_BA
    RHS[N+2]    -= gasup.enthalpy * mdot_BA * Vex / V_dV[1]
    Yup          = gasup.Y[1:end-1]
    RHS[3:N]    .+= mdot_BA * Yup
    RHS[N+3:2*N] .-= mdot_BA * Yup * Vex / V_dV[1]

    du .= Mass \ RHS
    return nothing
end

function f_intake_exhaust!(du, u, p, t)
    params, mode = p
    Mass = Mass_full(u, params, t)
    RHS  = params.scratch.RHS
    fill!(RHS, 0.0)
    RHS_Common!(RHS, u, params, t)
    N = params.gasses[1].gas.Nspec + 1

    V_dV         = params.Vfunc(t)
    TPX_exterior = mode.TPX_exterior

    # Use gas A slot as scratch for exterior state.
    # TPX_exterior must carry the full Nspec-element X vector.
    gas_ex = params.gasses[1]
    ct.setTPX(gas_ex, TPX_exterior)
    Pex     = TPX_exterior[2]
    K       = V_dV[2] > 0 ? params.K_vals[4] : params.K_vals[3]
    mdot_in = K * (Pex - u[N+1])   # positive = flow in; backflow allowed
    gasup   = mdot_in < 0 ? params.gasses[2] : gas_ex

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
        sol = ODE.solve(prob, alg, p=p, reltol=1e-8, abstol=1e-10, saveat=tsave, callback=cb)
    else
        sol = ODE.solve(prob, alg, p=p, reltol=1e-8, abstol=1e-10, saveat=tsave)
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
