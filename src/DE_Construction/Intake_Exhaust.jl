##########################################################
#=
Discrete intake/exhaust valve events.
The external volume always opens/closes when pressures equalize.

Exhaust: interior P ≥ exterior P → pressure drops to exterior, composition unchanged.
Intake:  interior P < exterior P → solve weighted mixture:
    P = P_intake
    ρ_tot*(h_new - h_cyl) = (Δm/V)*(h_ext - h_new)
    ρ_tot*(Y_new - Y_cyl) = (Δm/V)*(Y_ext - Y_new)
  where ρ_tot = ρ_cyl + Δm/V, and knowing P + ρ_tot + Y → h via Cantera.
=#
##########################################################

function open_valve_discrete!(u, TPX_exterior, gas, cfg::ModelConfig)
    if u[1] >= TPX_exterior[2]
        # Exhaust: pressure equalises, composition unchanged
        u[1] = TPX_exterior[2]
    else
        # Intake: solve for mixed state
        TPX_interior = (u[2], u[1], u[3:end])
        ct.setTPX(gas, TPX_interior, cfg.spec_ind)
        DHY_int = [gas.rho; gas.enthalpy; gas.Y[cfg.spec_ind[1:end-1]]]

        ct.setTPX(gas, TPX_exterior, cfg.spec_ind)
        DHY_ex = [gas.rho; gas.enthalpy; gas.Y[cfg.spec_ind[1:end-1]]]

        f!(du, x) = _valve_fsolve(du, x, (DHY_int, DHY_ex, gas, TPX_exterior[2], cfg))
        x0  = 0.5 * (DHY_int + DHY_ex)
        sol = NLsolve.nlsolve(f!, x0, method=:newton)

        ct.setHPY(gas, (sol.zero[2], TPX_exterior[2], sol.zero[3:end]), cfg.spec_ind)
        u .= [gas.P; gas.T; gas.X[cfg.spec_ind[1:end-1]]]
    end
end

function _valve_fsolve(F, u, param)
    # u = [Δm/V, h, Y...]
    # param = ([ρ, h, Y...]_cyl, [ρ, h, Y...]_ext, gas, P_target, cfg)
    cfg = param[5]
    ρt  = param[1][1] + u[1]
    N   = length(u)
    F[1] = param[1][1] * (u[2] - param[1][2]) - u[1] * (param[2][2] - u[2])
    for i in 1:N-2
        F[i+1] = param[1][1] * (u[i+2] - param[1][2+i]) - u[1] * (param[2][2+i] - u[i+2])
    end
    ct.setHPY(param[3], (u[2], param[4], u[3:end]), cfg.spec_ind)
    F[end] = ρt - param[3].rho
end
