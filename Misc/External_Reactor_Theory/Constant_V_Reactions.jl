topdir=abspath((@__DIR__)*"./../../")
using Pkg
Pkg.activate(topdir*"Kinetic_Model/.")

import OrdinaryDiffEq  as ODE
# need to build a constant volume solver (uck....)
# Wait! I already have one!!

include(topdir*"Kinetic_Model/External_Volume_3Region/DE_Construction/DE_Model.jl")
fluid_props=DE_Model.fluid_props
gas=DE_Model.initialize_ideal_gas("gri30.yaml")
specs=["CO", "CO2", "H2", "CH3OH", "H2O", "N2"]
DE_Model.set_gas_constants(gas,specs)
gas=DE_Model.initialize_ideal_gas("gri30.yaml")
param=(gas, 1000.0)
spec_ind=DE_Model.spec_ind

function const_V_reactor!(du, u, p, t)
    gas=p[1]
    M=DE_Model.Mass_A(gas)
    rhocat=p[2]
    XA=[u[3:end]; 1-sum(u[3:end])]
    # assume that catalyst is at wall temperature
    r=[0.0, 0.0, 0.0]

    r.= rhocat*[fluid_props.r1(u[2],u[1]/100e3,XA),
        fluid_props.r2(u[2],u[1]/100e3,XA), 
        fluid_props.r3(u[2],u[1]/100e3,XA)] # mol/s

    # first block is region A -> just reactions
    RHS = [0;0;(DE_Model.spec_MW[1:end-1]*1e-3).*[-r[1]+r[2];# production of CO in kg/s
            -r[2]-r[3];# production of CO2
            -2*r[1]-r[2]-3*r[3];# production of H2
            r[1]+r[3];# production of CH3OH
            r[2]+r[3]]] # production of H2O,
    return du .= M\RHS
end



P0=10.0:10:200
T0=220:10:300
comp="CO2:0.25, H2:0.75";
#comp="CO:0.3, H2:0.66, CO2:0.04";
Pf=Matrix{Float64}(undef, length(P0), length(T0))
Tf=similar(Pf)

using Printf
for i in eachindex(P0), j in eachindex(T0)
    @printf("Beginning P0=%.2fbar, T0=%.2fC...",P0[i], T0[j])
    gas.TPX=(T0[j]+273.15, P0[i]*100e3, comp)
    u0=[gas.P; gas.T; fluid_props.X(gas)[spec_ind[1:end-1]]]
    # Run kinetics through to comple
    prob=ODE.ODEProblem(const_V_reactor!,u0,(0,10))
    sol = ODE.solve(prob,ODE.RadauIIA5(autodiff=false),
        p=param, reltol=1e-8, abstol=1e-10)
    y=permutedims(stack(sol.u))
    Pf[i,j]=y[end,1]/100e3
    Tf[i,j]=y[end,2]
    @printf("Complete! Final P=%.2fbar\n",Pf[i,j])
end

begin
plot(P0, Pf, label=permutedims([@sprintf("T0=%.0fdegC",t) for t in T0]),
    xlabel="Initial P", ylabel="Final P")
plot!(P0,P0, linestyle=:dash, color=:black, label=:none)
end

begin
plot(T0, Tf'[:,2:2:end] .- 273.15, label=permutedims([@sprintf("P0=%.0fbar",t) for t in P0[2:2:end]]),
    xlabel="Initial T", ylabel="Final T", legend=:topleft)
plot!(T0,T0, linestyle=:dash, color=:black, label=:none)
end

begin
plot(sol.t,y[:,3:end], label=permutedims(specs[1:end-1]),
    xlabel="Time (s)", ylabel="Mole Fraction", legend=:right)
plot!([NaN],[NaN], color=:black, label="Pressure")
plot!(twinx(), sol.t, [y[:,1]/100e3,0*sol.t .+ P0[i]], legend=:none, ylabel="Pressure",
    color=[:black :grey], linestyle=[:solid :dash])
end

plot(sol.t, y[:,2])