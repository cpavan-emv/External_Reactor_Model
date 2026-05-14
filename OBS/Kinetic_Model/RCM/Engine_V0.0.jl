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
    V=(CR-1)/2*(cos(pi*t/tcomp)+1)+1e-8
    dV=(CR-1)/2*(-1)*pi/tcomp*sin(pi*t/tcomp)
    return V, dV
end

# setup the initial conditions and the parameters
CR=12
P0=10
tcomp=50e-3
Vfunc(t)=V_dV(t, CR, tcomp)
fluid_props.setTPX(gas1, (300.0, 100e3, "H2:0.75, CO2:0.25"))
fluid_props.setTPX(gas2, (300.0, 100e3, "H2:0.75, CO2:0.25"))
param=(gas1, gas2, Vfunc, rhocat, 250+273.15)

spec_ind=DE_construct.spec_ind
u0=vcat([P0*100e3, 300.0, 300.0], 
    fluid_props.X(gas1)[spec_ind[1:end-1]], 
    fluid_props.X(gas2)[spec_ind[1:end-1]])
du=copy(u0)

using OrdinaryDiffEq
begin
tmax=tcomp*1.9999#1.0#2*tcomp;
tp0=range(0.0,tmax,1001) # times for saving solution
tp=copy(tp0) # times for saving solution

# full implicit solution
# Radau seems to be the best choice
prob=ODEProblem(DE_construct.f!,u0,(0.0,tmax))
sol = solve(prob,RadauIIA5(autodiff=false),p=param, reltol=1e-8, abstol=1e-10, saveat=tp0)#
y=permutedims(stack(sol.u))

    for i in range(1,10)
        prob2=ODEProblem(DE_construct.f!,y[end,:],(tmax+(2*(i-1))*tcomp, tmax+2*i*tcomp))
        sol2 = solve(prob2,RadauIIA5(autodiff=false),p=param, reltol=1e-8, abstol=1e-10, saveat=2*i*tcomp .+ tp0)#
        y=vcat(y,permutedims(stack(sol2.u)))
        tp=vcat(tp, 2*i*tcomp .+ tp0)
    end

end



# IMEX is less efficient because the Mass matrix needs to be inverted regardless
# If I could get of the mass matrix inversion, then this would probably be better

# # split IMEX solution
# # may be more efficient than fully implicity
# prob_split=SplitODEProblem(DE_construct.f1!,DE_construct.f2!,u0,(0.0,tmax))
# sol_split = solve(prob_split,KenCarp4(autodiff=false),p=param, reltol=1e-8, abstol=1e-10, dtmax=5e-2,saveat=tp)#
# y_split=permutedims(stack(sol_split.u))


# using BenchmarkTools
# @btime solve(prob,RadauIIA5(autodiff=false),p=param, reltol=1e-8, abstol=1e-10, saveat=tp)#
# @btime solve(prob_split,KenCarp4(autodiff=false),p=param, reltol=1e-8, abstol=1e-10, dtmax=5e-2,saveat=tp)#
# @btime solve(prob_split,RadauIIA5(autodiff=false),p=param, reltol=1e-8, abstol=1e-10,saveat=tp)#

V=permutedims(stack(Vfunc.(tp)))[:,1]
dV=permutedims(stack(Vfunc.(tp)))[:,2]

fluid_props.setTPX(gas1, (300.0,100e3,gas1.X))
gamma=fluid_props.cp(gas1)/fluid_props.cv(gas1)

Nspec=length(specs)
#y=y_rxn
#y=y_comp
using Plots, Measures, LaTeXStrings
begin
p=[Plots.plot() for _ in 1:4]
xrng=(10*tcomp,12*tcomp*1e3)
p[1]=plot(tp*1e3,y[:,1]/100e3, 
    xlabel="Time (ms)", ylabel="Pressure (bar)", label="Pressure",
    leftmargin=10mm, title="Pressure", legend=:bottomright, xlim=xrng)

p[2]=plot(tp*1e3,y[:,2:3] .- 273.15, legend=:bottomright, label=["Region A" "Region B"],
    xlabel="Time (ms)", ylabel="Temp (deg C)",
    leftmargin=10mm, title="Temperature",xlim=xrng)

p[3]=plot(tp*1e3,y[:,4:4+4], label=permutedims(specs),
    xlabel="Time (ms)", ylabel="Mole Fraction", legend=:topleft,
    leftmargin=10mm, title="Composition, Region A",xlim=xrng)
plot!(p[3],tp*1e3, 1.0 .- sum(y[:,4:4+4], dims=2), label=specs[end],xlim=xrng)

p[4]=plot(tp*1e3,y[:,3+Nspec:end], label=permutedims(specs),
    xlabel="Time (ms)", ylabel="Mole Fraction", legend=:topleft,
    leftmargin=10mm, title="Composition, Region B",xlim=xrng)
plot!(p[4],tp*1e3, 1.0 .- sum(y[:,3+Nspec:end], dims=2), label=specs[end], xlim=xrng)


fig=plot(p..., size=(1080,720))
end
savefig(fig,topdir*"output/figures/Engine1")

#####################################################################
# Conservation checks
yeval=y
rho_uA=get_rhoU(gas1,yeval,1)
rho_uB=get_rhoU(gas2,yeval,2)
M=Mtot(rho_uA[:,1], rho_uB[:,1], V)
plot(tp,M)

U=Utot(rho_uA,rho_uB,V)
Ein=Ecomp(yeval[:,1],V)
plot(tp,U-Ein,legend=:bottomright)