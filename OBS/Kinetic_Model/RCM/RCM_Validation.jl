# this version moves all DE tools to separate 

# first get the environment set up
topdir=abspath((@__DIR__)*"/../")
using Pkg
Pkg.activate(topdir*"/src/.")
using LinearAlgebra

include(topdir*"src/fluid_properties/fluid_props.jl")
include(topdir*"src/DE_Construction/DE_Const.jl")
include(topdir*"src/evaluation/conservation_checks.jl")
DE_construct.fluid_props=fluid_props

# initialize the gas objects and user-set parameters
gas1=fluid_props.initialize_ideal_gas("gri30.yaml")
gas2=fluid_props.initialize_ideal_gas("gri30.yaml")
fluid_props.kingas=fluid_props.initialize_ideal_gas("gri30.yaml")
rhocat=1e3 # kg/m^3

specs=["CO", "CO2", "H2", "CH3OH", "H2O", "N2"]
DE_construct.set_gas_constants(gas1,specs)

# Define the compression function

function V_dV(t, CR, tcomp)
    # assume that V follows a cosine curve
    V=(CR-1)/2*(cos(pi*t/tcomp)+1)+1e-3
    dV=(CR-1)/2*(-1)*pi/tcomp*sin(pi*t/tcomp)
    if t>tcomp
        V=1e-3
        dV=0.0
    end
    return V, dV

end

# setup the initial conditions and the parameters
CR=8
P0=1
tcomp=50e-3
Vfunc(t)=V_dV(t, CR, tcomp)
gas_props0=(300.0, 100e3, "H2:0.75, CO2:0.25")
#gas_props0=(300.0, 100e3, "N2:1.0")
fluid_props.setTPX(gas1, gas_props0)
fluid_props.setTPX(gas2, gas_props0)
param=(gas1, gas2, Vfunc, rhocat)

spec_ind=DE_construct.spec_ind
u0=vcat([P0*100e3, 300.0, 300.0], 
    fluid_props.X(gas1)[spec_ind[1:end-1]], 
    fluid_props.X(gas2)[spec_ind[1:end-1]])
du=copy(u0)

using OrdinaryDiffEq
fluid_props.setTPX(gas1, (300.0,P0*100e3,gas1.X))
gamma1=fluid_props.cp(gas1)/fluid_props.cv(gas1)
fluid_props.setTPX(gas1, (550.0,P0*100e3,gas1.X))
gamma2=fluid_props.cp(gas1)/fluid_props.cv(gas1)
gamma_bar=0.5*(gamma1+gamma2)

########################################################################################
# Compression only validation - validates implementation of fluid model

tmax=4*tcomp;
tp=range(0.0,tmax,1001) # times for saving solution
V=permutedims(stack(Vfunc.(tp)))[:,1]
dV=permutedims(stack(Vfunc.(tp)))[:,2]

prob=ODEProblem(DE_construct.f2!,u0,(0.0,tmax))
sol = solve(prob,RadauIIA5(autodiff=false),p=param, reltol=1e-8, abstol=1e-10, saveat=tp)#
y=permutedims(stack(sol.u))


rho_uA=get_rhoU(gas1,y,1)
rho_uB=get_rhoU(gas2,y,2)
M=Mtot(rho_uA[:,1], rho_uB[:,1], V)
ΔM=M .- M[1]
U=Utot(rho_uA,rho_uB,V)
Ein=Ecomp(y[:,1],V)

using Plots, Measures, LaTeXStrings
begin
p=[Plots.plot() for _ in 1:3]
p[1]=plot(tp*1e3,u0[1]/100e3*( (1 .+ V[1])./(1 .+ V)).^gamma_bar, 
    xlabel="Time (ms)", ylabel="Pressure (bar)", label="Isentropic, fixed γ pressure", linewidth=2)
plot!(p[1],tp*1e3,y[:,1]/100e3, label="Detailed Simulation pressure", linestyle=:dash, linewidth=2)
plot!(p[1],[NaN], [NaN], color=:black, label="Temperature", leftmargin=10mm)

plot!(twinx(p[1]),tp*1e3,y[:,2] .- 273.15, legend=:none,
   ylabel="Temp (deg C)", color=:black, rightmargin=10mm)

p[2]=plot(tp*1e3,ΔM, legend=:none,
    xlabel="Time (ms)", ylabel=L"Change in Total Mass $(kg/V^A)$", leftmargin=10mm)

p[3]=plot(tp*1e3,[U Ein U .- Ein]/1e3, label=["Gas Energy" "Input Energy" "Gas - Input"],legend=:topright,
    xlabel="Time (ms)", ylabel=L"Total Energy $(kJ/V^A)$", leftmargin=10mm)

fig=plot(p..., layout=grid(3,1),size=(720,720))
end
savefig(fig,topdir*"output/figures/RCM_compression_validation")

########################################################################################
# Reaction only validation - validates implementation of reaction model
P0=100; T0=250+273.15
fluid_props.setTPX(gas1, gas_props0)
fluid_props.setTPX(gas2, gas_props0)
u0=vcat([P0*100e3, T0, T0], 
    fluid_props.X(gas1)[spec_ind[1:end-1]], 
    fluid_props.X(gas2)[spec_ind[1:end-1]])
Vfunc2(t)=V_dV(t .+ 2*tcomp, CR, tcomp)
param=(gas1, gas2, Vfunc2, rhocat)

tmax2=50

tp2=range(0.0,tmax2,1001) # times for saving solution
V2=permutedims(stack(Vfunc2.(tp2)))[:,1]


prob=ODEProblem(DE_construct.f1!,u0,(0.0,tmax2))
sol = solve(prob,RadauIIA5(autodiff=false),p=param, reltol=1e-8, abstol=1e-10, saveat=tp2)#
y=permutedims(stack(sol.u))


rho_uA=get_rhoU(gas1,y,1)
rho_uB=get_rhoU(gas2,y,2)
M=Mtot(rho_uA[:,1], rho_uB[:,1], V2)
ΔM=(M .- M[1]) ./ M[1]
U=Utot(rho_uA,rho_uB,V2)
ΔU=(U .- U[1]) ./ U[1]

using Plots; Measures, LaTeXStrings
begin
p=[Plots.plot() for _ in 1:4]
p[1]=plot(tp2,y[:,1]/100e3, label="Pressure",
    xlabel="Time (s)", ylabel="Pressure (bar)", linewidth=2)
plot!(p[1],[NaN], [NaN], color=:black, label="Temperature", leftmargin=10mm)

plot!(twinx(p[1]),tp2,y[:,2] .- 273.15, legend=:none,
   ylabel="Temp (deg C)", color=:black, rightmargin=10mm)

p[3]=plot(tp,y[:,4:4+4], label=permutedims(specs),
    xlabel="Time (s)", ylabel="Mole Fraction", legend=:topleft)
plot!(p[3],tp, 1.0 .- sum(y[:,4:4+4], dims=2), label=specs[end])

p[2]=plot(tp2,ΔM*100, legend=:none,
    xlabel="Time (s)", ylabel=L"Change in Total Mass $(\%)$", leftmargin=10mm)

p[4]=plot(tp2,ΔU*100,legend=:topright,
    xlabel="Time (s)", ylabel=L"Change in Internal Energy $(\%)$", leftmargin=10mm)



fig=plot(p..., layout=grid(2,2),size=(1080,720))
end
savefig(fig,topdir*"output/figures/RCM_rxn_validation")

########################################################################################
# Reacting

P0=10; T0=20+273.15
fluid_props.setTPX(gas1, gas_props0)
fluid_props.setTPX(gas2, gas_props0)
u0=vcat([P0*100e3, T0, T0], 
    fluid_props.X(gas1)[spec_ind[1:end-1]], 
    fluid_props.X(gas2)[spec_ind[1:end-1]])
param=(gas1, gas2, Vfunc, rhocat, 250+273.15)

tmax=20*tcomp;
tp=range(0.0,tmax,1001) # times for saving solution
V=permutedims(stack(Vfunc.(tp)))[:,1]
dV=permutedims(stack(Vfunc.(tp)))[:,2]

prob=ODEProblem(DE_construct.f!,u0,(0.0,tmax))
sol = solve(prob,RadauIIA5(autodiff=false),p=param, reltol=1e-8, abstol=1e-10, saveat=tp)#
y=permutedims(stack(sol.u))


rho_uA=get_rhoU(gas1,y,1)
rho_uB=get_rhoU(gas2,y,2)
M=Mtot(rho_uA[:,1], rho_uB[:,1], V)
ΔM=M .- M[1]
U=Utot(rho_uA,rho_uB,V)
Ein=Ecomp(y[:,1],V)

using Plots; Measures, LaTeXStrings
begin
p=[Plots.plot() for _ in 1:4]
p[1]=plot(tp*1e3,y[:,1]/100e3, label="Pressure",
    xlabel="Time (ms)", ylabel="Pressure (bar)", linewidth=2)
plot!(p[1],[NaN], [NaN], color=:black, label="Temperature", leftmargin=10mm)

plot!(twinx(p[1]),tp*1e3,y[:,2] .- 273.15, legend=:none,
   ylabel="Temp (deg C)", color=:black, rightmargin=10mm)

p[3]=plot(tp*1e3,y[:,4:4+4], label=permutedims(specs),
   xlabel="Time (ms)", ylabel="Mole Fraction", legend=:topleft)
plot!(p[3],tp, 1.0 .- sum(y[:,4:4+4], dims=2), label=specs[end])

p[2]=plot(tp*1e3,ΔM*100, legend=:none,
   xlabel="Time (s)", ylabel=L"Change in Total Mass $(\%)$", leftmargin=10mm)

p[4]=plot(tp*1e3,[U Ein U .- Ein]/1e3, label=["Gas Energy" "Input Energy" "Gas - Input"],legend=:topright,
   xlabel="Time (ms)", ylabel=L"Total Energy $(kJ/V^A)$", leftmargin=10mm)

fig=plot(p..., layout=grid(2,2),size=(1080,720))
end
savefig(fig,topdir*"output/figures/RCM_Simulation_2")


########################################################################################
# Approx. N2 compression

P0=1; T0=20+273.15
gas_props0=(300.0, 100e3, "H2:0.001, CO2:0.001, N2:0.998")
fluid_props.setTPX(gas1, gas_props0)
fluid_props.setTPX(gas2, gas_props0)
u0=vcat([P0*100e3, T0, T0], 
    fluid_props.X(gas1)[spec_ind[1:end-1]], 
    fluid_props.X(gas2)[spec_ind[1:end-1]])
param=(gas1, gas2, Vfunc, rhocat, 250+273.15)

tmax=20*tcomp;
tp=range(0.0,tmax,1001) # times for saving solution
V=permutedims(stack(Vfunc.(tp)))[:,1]
dV=permutedims(stack(Vfunc.(tp)))[:,2]

prob=ODEProblem(DE_construct.f2_demo!,u0,(0.0,tmax))
sol = solve(prob,RadauIIA5(autodiff=false),p=param, reltol=1e-8, abstol=1e-10, saveat=tp)#
y=permutedims(stack(sol.u))


rho_uA=get_rhoU(gas1,y,1)
rho_uB=get_rhoU(gas2,y,2)
M=Mtot(rho_uA[:,1], rho_uB[:,1], V)
ΔM=M .- M[1]
U=Utot(rho_uA,rho_uB,V)
Ein=Ecomp(y[:,1],V)

using Plots; Measures, LaTeXStrings
begin
p=[Plots.plot() for _ in 1:1]
p[1]=plot(tp*1e3,y[:,1]/100e3, label="Pressure",
    xlabel="Time (ms)", ylabel="Pressure (bar)", linewidth=2)
plot!(p[1],[NaN], [NaN], color=:black, label="Temperature", leftmargin=10mm, ylim=(0,20))

plot!(twinx(p[1]),tp*1e3,y[:,2] .- 273.15, legend=:none,
   ylabel="Temp (deg C)", color=:black, rightmargin=10mm)

fig=plot(p...,size=(1080,720))
end
savefig(fig,topdir*"output/figures/RCM_Simulation_3")