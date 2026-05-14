# this formulation of the DE uses a fully implicit form of the differential equation
# it has a singular, but constant, mass matrix
# equations:
# d(UA)/dt = mdot * hB --> d(rhoA * eA)/dt = 1/VA *mdot*hB``
# d(UB)/dt = - mdot * hB --> d(rhoB * eB) = 1/VB*(rhoB*eB*dVB - mdot*hB)
# d(MB)/dt + d(MA) = 0 
# d(MA)/dt = F
# d(mdot)/dt - dFdt = 0
# vars: 6
# equations: 5
# algebraic equation for closure
# PA=PB

# u vector order: [MA, MB, rho*uA, rho*uB, F, mdot]

function DAE_nonreactive!(du, u, p, t)
    gasA=p[1] # for passing to fluid props
    gasB=p[2] # for passing to fluid props
    V, dV = p[3](t) # function of t returning V and dV

    # first convert solution vector to gas state variables
    # take Normalize all volumes such that VA=1
    rhoA=u[1]
    rhoB=u[2]/V
    eA=u[3]/u[1]
    eB=u[4]/rhoB
    
    # now set the current state
    fluid_props.setUVX(gasA, (eA, 1/rhoA, gasA.X))
    fluid_props.setUVX(gasB, (eB, 1/rhoB,  gasB.X))

    # construct the equation
    du[1]=u[6]*fluid_props.enthalpy(gasB)
    du[2]=1/V*(rhoB*eB*dV-u[6]*fluid_props.enthalpy(gasB))
    du[3]=0.0 # in mass matrix
    du[4]=u[5]
    du[5]=0.0 # in mass matrix
    du[6]=fluid_props.P(gasA)-fluid_props.P(gasB) # algebraic equation
    return nothing
end

const mass_DAE_nonreactive = [
    0 0 1.0 0 0 0; # A energy
    0 0 0 1 0 0; # B energy
    1 1 0 0 0 0; # mass conservation
    1 0 0 0 0 0; # dummy 2nd order equation
    0 0 0 0 1 -1; # mdot
    0 0 0 0 0 0; # degenerate row for DAE
]

function J_numerical(fun!, u0, m=1e-8, b=1e-10)
    J=zeros(length(u0), length(u0))
    R0=similar(u0)
    fun!(R0,u0)
    R1=similar(R0)
    for i in eachindex(u0)
        u1=copy(u0)
        dui= u1[i]*m + b
        u1[i] += dui
        fun!(R1,u1)
        J[:,i] .= (R1-R0) / dui
    end
    return J
end


function DAE_nonreactive2!(res, du, u, p, t)
    gasA=p[1] # for passing to fluid props
    gasB=p[2] # for passing to fluid props
    V, dV = p[3](t) # function of t returning V and dV

    # first convert solution vector to gas state variables
    # take Normalize all volumes such that VA=1
    rhoA=u[1]
    rhoB=u[2]/V
    eA=u[3]/u[1]
    eB=u[4]/rhoB
    
    # now set the current state
    fluid_props.setUVX(gasA, (eA, 1/rhoA, gasA.X))
    fluid_props.setUVX(gasB, (eB, 1/rhoB,  gasB.X))

    # construct the equation
    RHS=similar(u)
    RHS[1]=u[6]*fluid_props.enthalpy(gasB)
    RHS[2]=1/V*(rhoB*eB*dV-u[6]*fluid_props.enthalpy(gasB))
    RHS[3]=0.0 # in mass matrix
    RHS[4]=u[5]
    RHS[5]=0.0 # in mass matrix
    RHS[6]=fluid_props.P(gasA)-fluid_props.P(gasB) # algebraic equation

    res .= mass_DAE_nonreactive * du - RHS
    # eliminate equation 5/ variable F

    return nothing
end

function implicit_Euler_step(f!, u0, param, t, dt,M)
    RHS=similar(u0)
    function fun!(res, x)
        f!(RHS, x, param, t)
        res.=(M*x-dt*RHS)-M*u0
    end
    sol=NLsolve.nlsolve(fun!, u0)
    return sol.zero
end

function implicit_Euler(f!, u0, param, Nstep, dt,M)
    u=Matrix{Float64}(undef, length(u0), Nstep+1)
    u[:,1].=u0
    t=0.0
    for i=2:Nstep+1
        t=t+dt
        u[:,i] = implicit_Euler_step(f!, u[:,i-1], param, t, dt, M)
    end
    return u
end

