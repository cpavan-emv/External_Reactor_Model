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

function V_dV(t, CR, param)
    # assume that V follows a cosine curve
    tcomp=param[1]
    tdelay=param[2]
    V=(CR-1)/2*(cos(pi*(t-tdelay)/tcomp)+1)+1e-3
    dV=(CR-1)/2*(-1)*pi/tcomp*sin(pi*(t-tdelay)/tcomp)
    if t>(tcomp+tdelay)
        V=1e-3
        dV=0.0
    elseif t<tdelay
        V=CR-1
        dV=0.0
    end
    return V, dV

end

# setup the initial conditions and the parameters
CR=8
P0=6
tcomp=10e-3
tdelay=10.0
Vfunc(t)=V_dV(t, CR, (tcomp, tdelay))
gas_props0=(250+273.15, 100e3, "H2:0.75, CO2:0.25, N2:0.0")
#gas_props0=(300.0, 100e3, "N2:1.0")
fluid_props.setTPX(gas1, gas_props0)
fluid_props.setTPX(gas2, gas_props0)
param=(gas1, gas2, Vfunc, rhocat, 250+273.15)

spec_ind=DE_construct.spec_ind
u0=vcat([P0*100e3, 250+273.15, 250+273.15], 
    fluid_props.X(gas1)[spec_ind[1:end-1]], 
    fluid_props.X(gas2)[spec_ind[1:end-1]])
du=copy(u0)

using OrdinaryDiffEq
tmax=tdelay+tcomp+0.5;
tp=range(0.0,tmax,10001) # times for saving solution

begin
prob=ODEProblem(DE_construct.f!,u0,(0.0,tdelay))
sol = solve(prob,RadauIIA5(autodiff=false),p=param, reltol=1e-8, abstol=1e-10, saveat=tp)#
y=permutedims(stack(sol.u))
u1=sol.u[end]
prob=ODEProblem(DE_construct.f!,u1,(tdelay,tmax))
sol = solve(prob,RadauIIA5(autodiff=false),p=param, reltol=1e-8, abstol=1e-10, saveat=tp)#
y=[y; permutedims(stack(sol.u))]
end




V=permutedims(stack(Vfunc.(tp)))[:,1]
dV=permutedims(stack(Vfunc.(tp)))[:,2]

fluid_props.setTPX(gas1, (300.0,100e3,gas1.X))
gamma=fluid_props.cp(gas1)/fluid_props.cv(gas1)


rho_uA=get_rhoU(gas1,y,1)
rho_uB=get_rhoU(gas2,y,2)
M=Mtot(rho_uA[:,1], rho_uB[:,1], V)

Xnet=(V.*y[:,9:end] .+ y[:,4:8])./(V .+ 1) # weighted average of the two region




using Plots; 
tplt=tp .- (tdelay+tcomp)
begin
p=[Plots.plot() for _ in 1:4]
p[1]=plot(tplt*1e3,u0[1]/100e3*( (1 .+ V[1])./(1 .+ V)).^gamma, 
    xlabel="Time (ms)", ylabel="Pressure (bar)", label="Isentropic, fixed γ")
plot!(p[1],tplt*1e3,y[:,1]/100e3, label="Detailed Simulation")

p[2]=plot(tplt*1e3,y[:,2] .- 273.15, legend=:none,
    xlabel="Time (ms)", ylabel="Temp (deg C)",)

p[3]=plot(tplt*1e3,Xnet, label=permutedims(specs),
    xlabel="Time (ms)", ylabel="Mole Fraction", legend=:topleft)
plot!(p[3],tplt*1e3, 1.0 .- sum(Xnet, dims=2), label=specs[end])

p[4]=plot(tplt*1e3,Xnet*100*1e4, label=permutedims(specs),
    xlabel="Time (ms)", ylabel="Concentration (ppm)", legend=:topleft,
    ylim=(0,1e4))

    for plt in p
        plot!(plt, xlim=(max(-tdelay, -10)*1e3,(tmax-tdelay)*1e3))
    end

fig=plot(p..., size=(1080,720))
end
#savefig(fig,topdir*"output/figures/heat_loss_validation.png")


#DF=DataFrame(:t=>tp, :P=>y[:,1])
#CSV.write("heat_loss_validation_data_segmented.csv", DF)

# #####################################################################
# # Conservation checks
yeval=y
rho_uA=get_rhoU(gas1,yeval,1)
rho_uB=get_rhoU(gas2,yeval,2)
M=Mtot(rho_uA[:,1], rho_uB[:,1], V)
plot(tp,M)

U=Utot(rho_uA,rho_uB,V)
Ein=Ecomp(yeval[:,1],V)
plot(tp,U-Ein,legend=:bottomright)

ΔM=(M .- M[1]) ./ M[1]
ΔU=(U .- U[1]) ./ U[1]
# u=y
# r=Matrix{Float64}(undef, length(tp), 3)
# for i in axes(r,1)
#     r[i,:] .= [fluid_props.r1(u[i,2],u[i,1]/100e3, [u[i,4:8];0.0]);
#     fluid_props.r2(u[i,2],u[i,1]/100e3, [u[i,4:8];0.0]);
#     fluid_props.r3(u[i,2],u[i,1]/100e3, [u[i,4:8];0.0])];

# end
# plot(tp,r, label=["CO hyd" "RWGS" "CO2 hyd"], xlim=[0, 300],
#     xlabel="Time (s)", ylabel="Reaction Rate (mol/s/kgcat)", legend=:topright)


using Plots, Measures, LaTeXStrings
begin
p=[Plots.plot() for _ in 1:4]
p[1]=plot(tp*1e3,y[:,1]/100e3, label="Pressure",
    xlabel="Time (ms)", ylabel="Pressure (bar)", linewidth=2)
plot!(p[1],[NaN], [NaN], color=:black, label="Temperature", leftmargin=10mm)

plot!(twinx(p[1]),tp*1e3,y[:,2] .- 273.15, legend=:none,
   ylabel="Temp (deg C)", color=:black, rightmargin=10mm)

p[3]=plot(tp*1e3,y[:,4:4+4], label=permutedims(specs),
    xlabel="Time (ms)", ylabel="Mole Fraction", legend=:topright)
plot!(p[3],tp, 1.0 .- sum(y[:,4:4+4], dims=2), label=specs[end])

p[2]=plot(tp*1e3,ΔM*100, legend=:none,
    xlabel="Time (ms)", ylabel=L"Change in Total Mass $(\%)$", leftmargin=10mm)

p[4]=plot(tp*1e3,[U Ein U .- Ein]/1e3, label=["Gas Energy" "Input Energy" "Gas - Input"],legend=:topright,
    xlabel="Time (ms)", ylabel=L"Total Energy $(kJ/V^A)$", leftmargin=10mm)


fig=plot(p..., layout=grid(2,2),size=(1080,720))
end
savefig(fig,topdir*"output/figures/RCM_Simulation_5.png")