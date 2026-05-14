topdir=abspath((@__DIR__)*"/../")
using Pkg
Pkg.activate(topdir*"/src/.")

include(topdir*"src/fluid_properties/fluid_props.jl")

gas1=fluid_props.initialize_ideal_gas("gri30.yaml")
gas2=fluid_props.initialize_ideal_gas("gri30.yaml")

function Mmat(gasA, gasB, Vstar)
    return [fluid_props.int_nrg(gasA)-fluid_props.enthalpy(gasB) 0 fluid_props.rho(gasA) 0;
        fluid_props.enthalpy(gasB)  Vstar*fluid_props.int_nrg(gasB) 0 Vstar*fluid_props.rho(gasB);
        1   Vstar   0   0]
end

function Dmat(gasA, gasB)
    return [fluid_props.rho(gasA)/fluid_props.P(gasA) -fluid_props.rho(gasA)/fluid_props.T(gasA)    0;
        fluid_props.rho(gasB)/fluid_props.P(gasB)   0   -fluid_props.rho(gasB)/fluid_props.T(gasB);
        0   fluid_props.cv(gasA)    0;
        0   0   fluid_props.cv(gasB)]
end

function V_dV(t, CompRatio, tcomp)
    if t<tcomp
        dV=-(CompRatio-1)/tcomp
        return (CompRatio-1) + dV*t, dV
    else
        return 0.0, 0.0
    end 
end

function f!(du, u, p, t)
    gasA=p[1] # for passing to fluid props
    gasB=p[2] # for passing to fluid props
    V, dV = p[3](t) # function of t returning V and dV
    fluid_props.setTPX(gasA, (u[2], u[1], gasA.X))
    fluid_props.setTPX(gasB, (u[3], u[1], gasB.X))
    M=Mmat(gasA, gasB, V)
    # M[:,2]./=-V
    D=Dmat(gasA, gasB)
    #D[:,3]./=-V ### <<<<<<<<<<WHHHHHHHHHHHYYYYYYYYYYY
    M2=(M*D)
    #M2^(-1)
    du .= M2\[0; -fluid_props.rho(gasB)*fluid_props.enthalpy(gasB)*dV; -fluid_props.rho(gasB)*dV]
    #du[end]*=-V
    return nothing
end

CR=16
Vfunc(t)=V_dV(t, CR, 10e-3)
fluid_props.setTPX(gas1, (300.0, 100e3, "N2:0.79, O2:0.21"))
fluid_props.setTPX(gas2, (300.0, 100e3, "N2:0.79, O2:0.21"))
param=(gas1, gas2, Vfunc)

u0=[100e3, 300.0, 300.0]
du=copy(u0)

f!(du, u0, param, 0.0)
using OrdinaryDiffEq
tmax=9.9e-3;
prob=ODEProblem(f!,u0,(0.0,tmax))
sol = solve(prob,Tsit5(),p=param)

tp=range(0.0,tmax,101)
y=permutedims(stack(sol.(tp)))
V=permutedims(stack(Vfunc.(tp)))[:,1]

fluid_props.setTPX(gas1, (300.0,100e3,gas1.X))
gamma=fluid_props.cp(gas1)/fluid_props.cv(gas1)

using Plots; 
begin
plot(tp,y[:,1])
plot!(tp,100e3*( (1 .+ V[1])./(1 .+ V)).^gamma)
plot!(tp,100e3*( y[:,2]./300.0).^(gamma/(gamma-1)))
plot!(tp,100e3*( y[:,3]./300.0).^(gamma/(gamma-1)), yscale=:log)
end

x=(1 .+ V[1])./(1 .+ V)
begin
    plot(x, 100e3*x.^gamma, xscale=:log, yscale=:log)
    plot!(x, y[:,1])
    plot!(x,100e3*( y[:,2]./300.0).^(gamma/(gamma-1)))
    plot!(x,100e3*( y[:,3]./300.0).^(gamma/(gamma-1)))
end



###########################################################
# Try implicit Euler with original equations - 
# verify change of vars, may also be a better formulation if I can find
# a DAE solver for non-constant mass matrix

function Mmat2(gasA, gasB, Vstar)
    return [-fluid_props.enthalpy(gasB) 0 1 0;
        fluid_props.enthalpy(gasB)  0 0 Vstar;
        1   Vstar   0   0;
        0 0 0 0]
end

function f2!(du,u, dt, uold, p, t)
    gasA=p[1] # for passing to fluid props
    gasB=p[2] # for passing to fluid props
    V, dV = p[3](t) # function of t returning V and dV
    fluid_props.setUVX(gasA, (u[3]/u[1], 1/u[1], gasA.X))
    fluid_props.setUVX(gasB, (u[4]/u[2], 1/u[2], gasB.X))
    M=Mmat2(gasA, gasB, V)
    
    du .= (M*u-dt*[0; 
        -fluid_props.rho(gasB)*fluid_props.enthalpy(gasB)*dV; 
        -fluid_props.rho(gasB)*dV;
        fluid_props.P(gasA)-fluid_props.P(gasB)])-M*uold

    return nothing
end


using NLsolve

ts=range(0.0,tmax, 101)
dt=ts[2]-ts[1]

fluid_props.setTPX(gas1, (300.0,100e3,gas1.X))
fluid_props.setTPX(gas2, (300.0,100e3,gas2.X))

u0=[fluid_props.rho(gas1), fluid_props.rho(gas2), 
    fluid_props.rho(gas1)*fluid_props.int_nrg(gas1), 
    fluid_props.rho(gas2)*fluid_props.int_nrg(gas2)]

uold=copy(u0)
u=copy(u0)

usol=zeros(length(ts),4)
usol[1,:]=u0

iend=101
for i in range(2, iend)
    u0=copy(usol[i-1,:])

    fun(F,u)=f2!(F, u, dt, usol[i-1,:], param, ts[i])
    sol=nlsolve(fun, u0)
    usol[i,:].=sol.zero
end

fluid_props.setUVX(gas1, (usol[iend,3]/usol[iend,1], 1/usol[iend,1], gas1.X))
fluid_props.setUVX(gas2, (usol[iend,4]/usol[iend,2], 1/usol[iend,2], gas2.X))


function f3!(res,du,u, p, t)
    gasA=p[1] # for passing to fluid props
    gasB=p[2] # for passing to fluid props
    V, dV = p[3](t) # function of t returning V and dV
    fluid_props.setUVX(gasA, (u[3]/u[1], 1/u[1], gasA.X))
    fluid_props.setUVX(gasB, (u[4]/u[2], 1/u[2], gasB.X))
    M=Mmat2(gasA, gasB, V)
    
    res .=M*du - [0; 
        -fluid_props.rho(gasB)*fluid_props.enthalpy(gasB)*dV; 
        -fluid_props.rho(gasB)*dV;
        fluid_props.P(gasA)-fluid_props.P(gasB)]

    return nothing
end

daeprob=DAEProblem(f3!, 0.1*ones(4), u0, (0.0,1e-10))#tmax))
solve(daeprob, DFBDF(autodiff=false),p=param)