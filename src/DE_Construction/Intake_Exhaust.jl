##########################################################
#=
Discrete intake/exhaust valve events.
The external volume always opens/closes when pressures equalize.

Exhaust: interior P >= exterior P -> pressure drops to exterior, composition unchanged.
Intake:  interior P < exterior P -> solve weighted mixture:
    P = P_intake
    rho_tot*(h_new - h_cyl) = (dm/V)*(h_ext - h_new)
    rho_tot*(Y_new - Y_cyl) = (dm/V)*(Y_ext - Y_new)
  where rho_tot = rho_cyl + dm/V, and knowing P + rho_tot + Y -> h via Cantera.
=#
##########################################################

function open_valve_discrete!(u, TPX_exterior, gas)
    if u[1] >= TPX_exterior[2]
        # Exhaust: pressure equalises, composition unchanged
        u[1] = TPX_exterior[2]
    else
        # Intake: solve for mixed state
        # u = [P, T, X_1..X_{Nspec-1}]; reconstruct full X with bath gas complement
        X_partial = u[3:end]
        X_full    = [X_partial; 1.0 - sum(X_partial)]
        ct.setTPX(gas, (u[2], u[1], X_full))
        DHY_int = [gas.rho; gas.enthalpy; gas.Y[1:end-1]]

        ct.setTPX(gas, TPX_exterior)
        DHY_ex = [gas.rho; gas.enthalpy; gas.Y[1:end-1]]

        f!(du, x) = _valve_fsolve(du, x, (DHY_int, DHY_ex, gas, TPX_exterior[2]))
        x0  = 0.5 * (DHY_int + DHY_ex)
        sol = NLsolve.nlsolve(f!, x0, method=:newton)

        Y_partial = sol.zero[3:end]
        Y_full    = [Y_partial; 1.0 - sum(Y_partial)]
        ct.setHPY(gas, (sol.zero[2], TPX_exterior[2], Y_full))
        u .= [gas.P; gas.T; gas.X[1:end-1]]
    end
end

function _valve_fsolve(F, u, param)
    # u = [dm/V, h, Y_1..Y_{Nspec-1}]
    # param = ([rho, h, Y...]_cyl, [rho, h, Y...]_ext, gas, P_target)
    rho_t = param[1][1] + u[1]
    N     = length(u)
    F[1]  = param[1][1] * (u[2] - param[1][2]) - u[1] * (param[2][2] - u[2])
    for i in 1:N-2
        F[i+1] = param[1][1] * (u[i+2] - param[1][2+i]) - u[1] * (param[2][2+i] - u[i+2])
    end
    Y_partial = u[3:end]
    Y_full    = [Y_partial; 1.0 - sum(Y_partial)]
    ct.setHPY(param[3], (u[2], param[4], Y_full))
    F[end] = rho_t - param[3].rho
end
