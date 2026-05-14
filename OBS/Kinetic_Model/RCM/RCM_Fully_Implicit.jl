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

# Define the complression function

function V_dV(t, CR, tcomp)
    # assume that V follows a cosine curve
    V=(CR-1)/2*(cos(pi*t/tcomp)+1)+1e-5
    dV=(CR-1)/2*(-1)*pi/tcomp*sin(pi*t/tcomp)
    if t>tcomp
        V=1e-5
        dV=0.0
    end
    return V, dV

end

# setup the initial conditions and the parameters
CR=12
P0=10
tcomp=50e-3
Vfunc(t)=V_dV(t, CR, tcomp)
fluid_props.setTPX(gas1, (300.0, 100e3, "H2:0.75, CO2:0.25"))
fluid_props.setTPX(gas2, (300.0, 100e3, "H2:0.75, CO2:0.25"))
param=(gas1, gas2, Vfunc, rhocat)

V0,_=Vfunc(0)

spec_ind=DE_construct.spec_ind
u0=[fluid_props.rho(gas1);
    fluid_props.rho(gas2)*V0;
    fluid_props.int_nrg(gas1);
    fluid_props.int_nrg(gas2);
    0.0;
    0.0]
u0[3]*=u0[1]
u0[4]*=fluid_props.rho(gas2)
du=copy(u0)

using OrdinaryDiffEq
f = ODEFunction(DE_construct.DAE_nonreactive!, mass_matrix = DE_construct.mass_DAE_nonreactive)

tmax=0.2*tcomp;
tp=range(0.0,tmax,1001) # times for saving solution

f!(du, u, p, t)=DE_construct.DAE_nonreactive!(du, u, p, t)
f!(du, u0, param, 0.0)

du0=0*u0
prob = DAEProblem(DE_construct.DAE_nonreactive2!, du0,u0,(0.0,tmax))
sol=solve(prob, DImplicitEuler(autodiff=false), p=param)

M=DE_construct.mass_DAE_nonreactive
using LinearAlgebra
tmp=diagm(ones(12))
tmp[1:6,1:6] .= M
b=vcat(zeros(6), u0)
# I know initial value for everything except F
# I know d/dt initial value for everything 


fun!(du, u)=DE_construct.DAE_nonreactive!(du,u, param, 0)
J=DE_construct.J_numerical(fun!, u0)
tmp[1:6, 7:end] .= J

# full implicit solution
# Radau seems to be the best choice
prob=ODEProblem(f,u0,(0.0,tmax), :ShampineCollocationInit)
#sol=solve(prob, RadauIIA5(autodiff=false), p=param, reltol=1e-8)
sol=solve(prob, Rosenbrock23(autodiff=false), p=param)#, reltol=1e-8, saveat=tp)
y=permutedims(stack(sol.u))

tloc=sol.t

plot(tloc,y)

sol = solve(prob,RadauIIA5(autodiff=false),p=param, reltol=1e-8, abstol=1e-10, saveat=tp)#
y=permutedims(stack(sol.u))

# IMEX is less efficient because the Mass matrix needs to be inverted regardless
# If I could get of the mass matrix inversion, then this would probably be better

# # split IMEX solution
# # may be more efficient than fully implicity
# prob_split=SplitODEProblem(DE_construct.f1!,DE_construct.f2!,u0,(0.0,tmax))
# sol_split = solve(prob_split,KenCarp4(autodiff=false),p=param, reltol=1e-8, abstol=1e-10, dtmax=5e-2,saveat=tp)#
# y_split=permutedims(stack(sol_split.u))


# using BenchmarkTools
@btime solve(prob,RadauIIA5(autodiff=false),p=param, reltol=1e-8, abstol=1e-10, saveat=tp)#
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

#p[2]=plot(tp*1e3,y[:,2:3] .- 273.15, label=["Rxn Volume" "Disp. Volume"],
#    xlabel="Time (ms)", ylabel="Temp (deg C)",)
p[2]=plot(tp*1e3,y[:,2] .- 273.15, legend=:none,
    xlabel="Time (ms)", ylabel="Temp (deg C)",)

p[3]=plot(tp*1e3,y[:,4:4+4], label=permutedims(specs),
    xlabel="Time (ms)", ylabel="Mole Fraction", legend=:topleft)
plot!(p[3],tp*1e3, 1.0 .- sum(y[:,4:4+4], dims=2), label=specs[end])

fig=plot(p..., size=(1080,720))
end
savefig(fig,topdir*"output/figures/rxn_prelim4")

#####################################################################
# Conservation checks
yeval=y
rho_uA=get_rhoU(gas1,yeval,1)
rho_uB=get_rhoU(gas2,yeval,2)
plot(tp,V.*rho_uB[:,1] + rho_uA[:,1])
#plot(tp,rho_uB[:,1] + rho_uA[:,1])

#y=yeval

V=permutedims(stack(Vfunc.(tp)))[:,1]
# energy Conservation
Ein=similar(tp)
Ein[1]=0.0
deltaE = 0.5*(yeval[2:end,1]+yeval[1:end-1,1]) .* (V[2:end]-V[1:end-1])
for i in eachindex(tp[1:end-1])
    Ein[i+1] =Ein[i]+deltaE[i]
end
Ein*=-1

begin
plot(tp,V.*rho_uB[:,2].*rho_uB[:,1] + rho_uA[:,2].*rho_uA[:,1]-Ein,legend=:bottomright)
#plot!(tp,Ein)
end