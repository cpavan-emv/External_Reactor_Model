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
tcomp=10e-3
Vfunc(t)=V_dV(t, CR, tcomp)
gas_props0=(300.0, 100e3, "H2:0.01, CO2:0.01, N2:0.98")
#gas_props0=(300.0, 100e3, "N2:1.0")
fluid_props.setTPX(gas1, gas_props0)
fluid_props.setTPX(gas2, gas_props0)
param=(gas1, gas2, Vfunc, rhocat, 250+273.15)

spec_ind=DE_construct.spec_ind
u0=vcat([P0*100e3, 300.0, 300.0], 
    fluid_props.X(gas1)[spec_ind[1:end-1]], 
    fluid_props.X(gas2)[spec_ind[1:end-1]])
du=copy(u0)

using OrdinaryDiffEq
tmax=1.0#10*tcomp;
tp=range(0.0,tmax,1001) # times for saving solution

# full implicit solution
# Radau seems to be the best choice
#prob=ODEProblem(DE_construct.f!,u0,(0.0,tmax))
# non-reactive only
prob=ODEProblem(DE_construct.f!,u0,(0.0,tmax))
sol = solve(prob,RadauIIA5(autodiff=false),p=param, reltol=1e-8, abstol=1e-10, saveat=tp)#
y=permutedims(stack(sol.u))

# IMEX is less efficient because the Mass matrix needs to be inverted regardless
# If I could get of the mass matrix inversion, then this would probably be better

# # split IMEX solution
# # may be more efficient than fully implicity
# prob_split=SplitODEProblem(DE_construct.f1!,DE_construct.f2!,u0,(0.0,tmax))
# sol_split = solve(prob_split,KenCarp4(autodiff=false),p=param, reltol=1e-8, abstol=1e-10, dtmax=5e-2,saveat=tp)#
# y_split=permutedims(stack(sol_split.u))


using BenchmarkTools
# @btime solve(prob,RadauIIA5(autodiff=false),p=param, reltol=1e-8, abstol=1e-10, saveat=tp)#
# @btime solve(prob_split,KenCarp4(autodiff=false),p=param, reltol=1e-8, abstol=1e-10, dtmax=5e-2,saveat=tp)#
# @btime solve(prob_split,RadauIIA5(autodiff=false),p=param, reltol=1e-8, abstol=1e-10,saveat=tp)#

V=permutedims(stack(Vfunc.(tp)))[:,1]
dV=permutedims(stack(Vfunc.(tp)))[:,2]

fluid_props.setTPX(gas1, (300.0,100e3,gas1.X))
gamma=fluid_props.cp(gas1)/fluid_props.cv(gas1)


#y=y_rxn
#y=y_comp
using Plots; 
begin
p=[Plots.plot() for _ in 1:3]
p[1]=plot(tp*1e3,u0[1]/100e3*( (1 .+ V[1])./(1 .+ V)).^gamma, 
    xlabel="Time (ms)", ylabel="Pressure (bar)", label="Isentropic, fixed γ")
plot!(p[1],tp*1e3,y[:,1]/100e3, label="Detailed Simulation")

p[2]=plot(tp*1e3,y[:,2] .- 273.15, legend=:none,
    xlabel="Time (ms)", ylabel="Temp (deg C)",)

p[3]=plot(tp*1e3,y[:,4:4+4], label=permutedims(specs),
    xlabel="Time (ms)", ylabel="Mole Fraction", legend=:topleft)
plot!(p[3],tp*1e3, 1.0 .- sum(y[:,4:4+4], dims=2), label=specs[end])

fig=plot(p..., size=(1080,720))
end
savefig(fig,topdir*"output/figures/heat_loss_validation.png")


DF=DataFrame(:t=>tp, :P=>y[:,1])
CSV.write("heat_loss_validation_data_segmented.csv", DF)

# #####################################################################
# # Conservation checks
# yeval=y
# rho_uA=get_rhoU(gas1,yeval,1)
# rho_uB=get_rhoU(gas2,yeval,2)
# M=Mtot(rho_uA[:,1], rho_uB[:,1], V)
# plot(tp,M)

# U=Utot(rho_uA,rho_uB,V)
# Ein=Ecomp(yeval[:,1],V)
# plot(tp,U-Ein,legend=:bottomright)