##########################################################
#=
This is for intake and exhaust
    (I assume the external volume always opens/closes when pressure is equalized)
Exhaust is easy - pressure instantly drops to exhaust P with no change in composition or temp
Intake is a bit tricky, must solve:
    - P=Pintake
    - htot*rho_tot=hcyl*rho_cyl + Δm * hin / Vcyl
    - Ytot*rho_tot = Ycyl*rho_cyl + Δm * Yin / Vcyl
Unknowns are htot, Δm, rhotot, Ytot
BUT: rho_tot = rho_cyl + Δm/Vcyl and knowing P, rhotot and Y determines h



=#
#########################################################

function open_valve_discrete!(u,TPX_exterior, gas)
    if u[1]>=TPX_exterior[2]
        # flow is out the valve
        # only thing that matters is exterior pressure
        u[1]=TPX_exterior[2]
    else
        # flow is into the chamber
        # need weighted average of existing gas and new gas
        TPX_interior=(u[2],u[1],u[3:end])
        fluid_props.setTPX(gas,TPX_interior, spec_ind)
        DHY_int=[fluid_props.rho(gas); fluid_props.enthalpy(gas); fluid_props.Y(gas)[spec_ind[1:end-1]]]
        fluid_props.setTPX(gas,TPX_exterior, spec_ind)
        DHY_ex=[fluid_props.rho(gas); fluid_props.enthalpy(gas); fluid_props.Y(gas)[spec_ind[1:end-1]]]
        f!(du,u) = _valve_fsolve(du,u,(DHY_int,DHY_ex, gas, TPX_exterior[2]))
        x0=0.5*(DHY_int+DHY_ex)
        sol=NLsolve.nlsolve(f!, x0, method=:newton)
        fluid_props.setHPY(gas, (sol.zero[2], TPX_exterior[2], sol.zero[3:end]), spec_ind)
        u .= [fluid_props.P(gas); fluid_props.T(gas); fluid_props.X(gas)[spec_ind[1:end-1]]]
    end
end

function _valve_fsolve(F,u,param)
    # u is a vector of [Δm/V, h, Y]
    # param is ([ρ, h, Y]_cyl, [ρ, h, Y]_ex, gas, Ptarget) 
    ρt=param[1][1]+u[1]
    N=length(u)
    # energy balance
    F[1] = param[1][1]*(u[2]-param[1][2]) - u[1]*(param[2][2]-u[2])
    # species mass balance
    for i in range(1,N-2)
        F[i+1] = param[1][1]*(u[i+2]-param[1][2+i])-u[1]*(param[2][2+i]-u[i+2])
    end
    # pressure matching
    fluid_props.setHPY(param[3], (u[2], param[4], u[3:end]), spec_ind)
    F[end] = ρt - fluid_props.rho(param[3])
end


