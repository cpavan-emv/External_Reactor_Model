# first get the environment set up
topdir=abspath((@__DIR__)*"/../")
using Pkg
Pkg.activate(topdir*".")
using LinearAlgebra

include(topdir*"External_Volume_3Region/DE_Construction/DE_Model.jl")
include((@__DIR__)*"/conservation_checks.jl")
fluid_props=DE_Model.fluid_props

gas=DE_Model.initialize_ideal_gas("gri30.yaml")
specs=["CO", "CO2", "H2", "CH3OH", "H2O", "N2"]
DE_Model.set_gas_constants(gas,specs)


gasses=[DE_Model.initialize_ideal_gas("gri30.yaml") for _ in 1:3]  
rhocat=1e3 # kg/m^3 catalyst effective density


# Define the compression function
# V is defined as the displaced volume divided by the clearance volume
# V will go to 0 when piston is at TDC

Vdisp=1; # piston displacement (reference dimension)
Vext=1; # external volume
CR0=18 # Compression ratio with valve closed
# CR=(Vdisp+Vclear)/Vclear
Vclear=Vdisp/(CR0-1)

function V_dV(t, CR, param)
    # assume that V follows a cosine curve
    tcomp=param[1]
    tdelay=param[2]
    V=(CR-1)/2*(cos(pi*(t-tdelay)/tcomp)+1)+1e-3
    dV=(CR-1)/2*(-1)*pi/tcomp*sin(pi*(t-tdelay)/tcomp)
    return V, dV
end

# setup the initial conditions and the parameters
Pin=6
Pout=1
Preact=10

tcomp=50e-3
tdelay=0
Vfunc(t)=V_dV(t, CR0, (tcomp, tdelay))
gas_props0=(250+273.15, 100e3, "H2:0.75, CO2:0.25, CO:0.0")
gas_props1=(300.0, 100e3, "N2:1.0")
[fluid_props.setTPX(g, gas_props1) for g in gasses[2:3]] # initialize everything to nitrogen

#fluid_props.setTPX(gasA, gas_props1)
#fluid_props.setTPX(gasB, gas_props0)
#param_DC=(gasA, gasB, Vfunc_decoup, rhocat, 250+273.15)
#param_CC=(gasA, gasB, Vfunc_coup, rhocat, 250+273.15)
param=(gasses, Vfunc, rhocat, (250+273.15)*ones(3), 50/1000*ones(3), [0.002,0.002],1)

spec_ind=DE_Model.spec_ind
u0_ext=[[Preact*100e3, 250+273.15];fluid_props.X(gasses[1])[spec_ind[1:end-1]]]
u0_cyl=[[Pin*100e3, 25+273.15];fluid_props.X(gasses[2])[spec_ind[1:end-1]]]
u0=[u0_ext;u0_cyl;u0_cyl]

gas_tmp=gasses[1]
fluid_props.setTPX(gas_tmp, gas_props0)

TPX_intake=(25+273.15, Pin*100e3,fluid_props.X(gas_tmp)[spec_ind[1:end-1]])
TPX_exhaust=(25+273.15, Pout*100e3,fluid_props.X(gas_tmp)[spec_ind[1:end-1]])
du=copy(u0)

param=(gasses, Vfunc, rhocat, [250.0,85.0,85.0] .+ 273.15, [50,500,500]*1e-3, [0.002,0.002, 0.02, 0.02],1)
param_in=(param[1:end-1]...,TPX_intake)
param_ex=(param[1:end-1]...,TPX_intake)


tmax=tcomp*2
topen=tcomp/10
tp=range(0.0,tmax,1001) # times for saving solution (each cycle)
tex=tp[1:floor(Int, length(tp)/2)+1]
tin=tp[floor(Int, length(tp)/2)+1:end]
Nspec=DE_Model.Nspec

##################
Ncycle=30
tcycle=tcomp*2

u=u0';
t=0.0;

for N in 1:Ncycle
    tstart=t[end]

# Starting at BDC (start of exhaust stroke)
tloc,y=DE_Model.evolve_intake_exhaust(u[end,:],param_ex,tex)
t=[t;tloc.+tstart]
u=[u;y]
tloc,y=DE_Model.evolve_intake_exhaust(u[end,:],param_in,tin)
t=[t;tloc.+tstart]
u=[u;y]
tstart=t[end]

# Next is compression - compress in isolation until the pressure in chambers equalizes
cond(u,t,int)=DE_Model.Vars_Eq_condition(u,t,int,[1,Nspec+2])
tloc,y=DE_Model.evolve_decoupled(u[end,:],param,tp,cond)
t=[t;tloc.+tstart]
u=[u;y]


# once the pressures equilizes, chambers are connected 
# they stay connected for a fixed time past TDC
tp2=[tloc[end];tp[tp .> tloc[end]]]
tp2=[tp2[tp2.<tcomp+topen];tcomp+topen]
tloc,y=DE_Model.evolve_coupled(u[end,:],param,tp2)
t=[t;tloc.+tstart]
u=[u;y]

# Expansion continues with chambers disconnected
tp3=[tp2[end];tp[tp .> tp2[end]]]
tloc,y=DE_Model.evolve_decoupled(u[end,:],param,tp3)
t=[t;tloc.+tstart]
u=[u;y]

# finally, the exhaust valve opens
#t=[t;t[end]]
#u=[u;u[end,:]']
#DE_Model.open_valve_discrete!((@view u[end,Nspec+2:end]), TPX_exhaust, gasB)
    println("Cycle $N Complete")
end


yA=u[:,1:1+Nspec]
yB=u[:,Nspec+2:2*(1+Nspec)]
yC=u[:,2*(1+Nspec)+1:end]

using Plots, Measures
begin
    plt=[Plots.plot() for _ in 1:4]
    plot!(plt[1],t, [yA[:,2] yB[:,2] yC[:,2]] .- 273, xlabel="Time (s)", ylabel="Temp (degC)",
        label=["Ext. Vol" "Clear. Vol" "Disp. Vol"])
    plot!(plt[2],t, [yA[:,1] yB[:,1] yC[:,1]]/100e3, xlabel="Time (s)", ylabel="Pressure (bar)",
        label=["Ext. Vol" "Clear. Vol" "Disp. Vol"])
    XA=yA[:,3:end]
    XA=[XA 1 .- sum(XA, dims=2)]
    XC=yC[:,3:end]
    XC=[XC 1 .- sum(XC, dims=2)]
    plot!(plt[3],t, XA*100, xlabel="Time (s)", ylabel="Mole Fraction (ext. vol)", ylim=(0,1),
        label=permutedims(specs))
    plot!(plt[4],t, XC*100, xlabel="Time (s)", ylabel="Mole Fraction (disp. vol)", ylim=(0,1),
        label=permutedims(specs))
    for p in plt
        plot!(p, xlim=(52*tcomp, 60*tcomp))
    end
    plot(plt..., size=(1080,720), margin=10mm)
end






prob=ODE.ODEProblem(DE_Model.f_CC!,u0,(0.0,tmax))
sol = ODE.solve(prob,RadauIIA5(autodiff=false),p=param, reltol=1e-8, abstol=1e-10, saveat=tp)#
y=permutedims(stack(sol.u))

# test event
function condition(u,t,integrator)
    return u[1] - 10*P0*100e3
end

cb=ODE.ContinuousCallback(condition, (integrator)->ODE.terminate!(integrator))

prob=ODE.ODEProblem(DE_Model.f_CC!,u0,(0.0,tmax))
sol = ODE.solve(prob,RadauIIA5(autodiff=false),p=param, reltol=1e-8, abstol=1e-10, saveat=tp, callback=cb)#
y=permutedims(stack(sol.u))

V=permutedims(stack(Vfunc.(sol.t)))[:,1]
dV=permutedims(stack(Vfunc.(sol.t)))[:,2]

fluid_props.setTPX(gasA, (300.0,100e3,fluid_props.X(gasA)))
gamma=fluid_props.cp(gasA)/fluid_props.cv(gasA)

rho_uA=get_rhoU(gasA,y,1)
rho_uB=get_rhoU(gasB,y,2)
M=Mtot(rho_uA[:,1], rho_uB[:,1], V)

Xnet=(V.*y[:,9:end] .+ y[:,4:8])./(V .+ 1) # weighted average of the two region

u0_CV=vcat([P0*10*100e3, 250+273.15], 
    fluid_props.X(gasA)[spec_ind[1:end-1]])
prob=ODE.ODEProblem(DE_Model.f_CV!,u0_CV,(0.0,tmax))
sol = ODE.solve(prob,RadauIIA5(autodiff=false),p=param, reltol=1e-8, abstol=1e-10, saveat=tp)#
y=permutedims(stack(sol.u))
rho_uA=get_rhoU_single_vol(gasA,y)


u0_Comp=vcat([P0*100e3, 250+273.15], 
    fluid_props.X(gasB)[spec_ind[1:end-1]])
prob=ODE.ODEProblem(DE_Model.f_Comp!,u0_Comp,(0.0,tmax))
sol = ODE.solve(prob,RadauIIA5(autodiff=false),p=param, reltol=1e-8, abstol=1e-10, saveat=tp)#
y=permutedims(stack(sol.u))
V=permutedims(stack(Vfunc.(sol.t)))[:,1]
rho_uB=get_rhoU_single_vol(gasB,y)[:,1] .* V

using Plots; 
tplt=sol.t#tp #.- (tdelay+tcomp)
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
        plot!(plt, xlim=(max(-tdelay, -10)*1e3,tplt[end]*1e3))
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