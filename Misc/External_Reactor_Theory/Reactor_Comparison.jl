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
rhocat_bulk=1000.0
param=(gas, rhocat_bulk)
spec_ind=DE_Model.spec_ind


Pr=100*100e3
Tr=290+273.15
comp="CO2:0.25, H2:0.75";

function const_V_reactor_adiabatic!(du, u, p, t)
    gas=p[1]
    fluid_props.setTPX(gas, (u[2], u[1], u[3:end]), spec_ind)
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

function const_V_reactor_isothermal!(du, u, p, t)
    gas=p[1]
    fluid_props.setTPX(gas, (u[2], u[1], u[3:end]), spec_ind)
    M=DE_Model.Mass_A(gas)
    rhocat=p[2]
    Tr=p[3]
    XA=[u[3:end]; 1-sum(u[3:end])]
    # assume that catalyst is at wall temperature
    r=[0.0, 0.0, 0.0]

    r.= rhocat*[fluid_props.r1(Tr,u[1]/100e3,XA),
        fluid_props.r2(Tr,u[1]/100e3,XA), 
        fluid_props.r3(Tr,u[1]/100e3,XA)] # mol/s

    # add massive heat loss to keep isothermal

    # first block is region A -> just reactions
    RHS = [0;-1e9*(u[2]-Tr);(DE_Model.spec_MW[1:end-1]*1e-3).*[-r[1]+r[2];# production of CO in kg/s
            -r[2]-r[3];# production of CO2
            -2*r[1]-r[2]-3*r[3];# production of H2
            r[1]+r[3];# production of CH3OH
            r[2]+r[3]]] # production of H2O,
    return du .= M\RHS
end


# Get STP conditions
gas.TPX=(273.15, 101.325e3, comp)
ρSTP=gas.density
tsave=0:0.1:100
GHSV_constV_i=Vector{Float64}(undef, length(tsave))
Xall_constV_i=Matrix{Float64}(undef, length(tsave), length(specs))
GHSV_constV_a=similar(GHSV_constV_i)
Xall_constV_a=similar(Xall_constV_i)


gas.TPX=(Tr, Pr, comp)
VSTP=gas.density/ρSTP
u0=[gas.P; gas.T; fluid_props.X(gas)[spec_ind[1:end-1]]]
param_constV=(gas,rhocat_bulk, Tr)

prob1=ODE.ODEProblem(const_V_reactor_adiabatic!,u0,(tsave[1], tsave[end]),saveat=tsave)
sol1 = ODE.solve(prob1,ODE.RadauIIA5(autodiff=false),
    p=param_constV, reltol=1e-8, abstol=1e-10)

prob2=ODE.ODEProblem(const_V_reactor_isothermal!,u0,(tsave[1], tsave[end]),saveat=tsave)
sol2 = ODE.solve(prob2,ODE.RadauIIA5(autodiff=false),
    p=param_constV, reltol=1e-8, abstol=1e-10)

y1=permutedims(stack(sol1.u))
y2=permutedims(stack(sol2.u))
GHSV_constV_a.=VSTP ./ (rhocat_bulk*sol1.t[:])*3600 # L/gcat/hr
Xall_constV_a .= [y1[:, 3:end] (1 .- sum(y1[:,3:end], dims=2))]
GHSV_constV_i.=VSTP ./ (rhocat_bulk*sol2.t[:])*3600 # L/gcat/hr
Xall_constV_i .= [y2[:, 3:end] (1 .- sum(y2[:,3:end], dims=2))]

using Plots, Measures, Printf, LaTeXStrings, ColorSchemes
cols=ColorSchemes.tab10
begin
plt_constV=plot(GHSV_constV_a*1e3,Xall_constV_a, label=permutedims(specs), 
    color=permutedims(cols[1:length(specs)]))
plot!(plt_constV,GHSV_constV_i*1e3,Xall_constV_i, label=:none, 
    color=permutedims(cols[1:length(specs)]),
    linestyle=:dot)
plot!(plt_constV,xaxis=:log, xlim=(1e3,1e6),ylabel="Mole Fraction",xlabel="GHSV (scc/gcat/hr)",
    title="Constant V Reactor\n Solid=Adiabatic, Dash=Isothermal", ylim=(0,0.30))
plot!(plt_constV,[NaN],[NaN], label="T (ad)", color=:black, linestyle=:solid)
plot!(twinx(plt_constV), GHSV_constV_a*1e3, y1[:,2] .- 273.15, ylabel="T (degC)",
    color=:black, linestyle=:solid, xlim=(1e3,1e6), legend=:none, xaxis=:log, ylim=(240, 300))
end

#################################################
# Now do PFR
# use formulation I alread developed (python model)
param_PFR=(gas, rhocat_bulk, ρSTP, GHSV_constV_a[end]/3600)
zstar=range(0, 1, length(GHSV_constV_a))

# isothermal
function PF_reactor_isothermal!(du, u, p, t)
    gas=p[1]
    fluid_props.setTPX(gas, (Tr, Pr, u[2:end]), spec_ind)
    M=DE_Model.Bcom(gas)[3:end, 3:end]
    rhocat=p[2]
    rho_STP=p[3]
    GHSV_end=p[4]
    XA=[u[2:end]; 1-sum(u[2:end])]
    # assume that catalyst is at wall temperature
    r=[0.0, 0.0, 0.0]

    r.= [fluid_props.r1(Tr,Pr/100e3,XA),
        fluid_props.r2(Tr,Pr/100e3,XA), 
        fluid_props.r3(Tr,Pr/100e3,XA)] # mol/s

    # first block is region A -> just reactions
    RHS = (DE_Model.spec_MW[1:end-1]*1e-3).*[-r[1]+r[2];# production of CO in kg/s
            -r[2]-r[3];# production of CO2
            -2*r[1]-r[2]-3*r[3];# production of H2
            r[1]+r[3];# production of CH3OH
            r[2]+r[3]] # production of H2O,
    du[2:end] .= M\(RHS./(rho_STP*GHSV_end))
    du[1] = gas.density # this is the equation for getting residence time
end

u0_PFR_i=[0;u0[3:end]]
prob3=ODE.ODEProblem(PF_reactor_isothermal!,u0_PFR_i,(0,1),saveat=zstar)
sol3 = ODE.solve(prob3,ODE.RadauIIA5(autodiff=false),
    p=param_PFR, reltol=1e-8, abstol=1e-10)
y_PFR_i=permutedims(stack(sol3.u))
GHSV_PFR_i=GHSV_constV_a[end]*y_PFR_i[end,1] ./ y_PFR_i[:,1]
Xall_PFR_i=[y_PFR_i[:,2:end] (1 .- sum(y_PFR_i[:,2:end], dims=2))]


# p_PFR=plot(GHSV_PFR*1e3,Xall_PFR, label=permutedims(specs), xlabel="GHSV (scc/gcat/hr)",ylabel="Mole Fraction",
#     xaxis=:log, xlim=(1e4,1e6))


# adiabatic PFR

function PF_reactor_adiabatic!(du, u, p, t)
    gas=p[1]
    fluid_props.setHPX(gas, (u[2], Pr, u[3:end]), spec_ind)
    M=DE_Model.Bcom(gas)[3:end, 3:end]
    T=gas.T
    rho_STP=p[3]
    GHSV_end=p[4]
    XA=[u[3:end]; 1-sum(u[3:end])]
    # assume that catalyst is at wall temperature
    r=[0.0, 0.0, 0.0]

    r.= [fluid_props.r1(T,Pr/100e3,XA),
        fluid_props.r2(T,Pr/100e3,XA), 
        fluid_props.r3(T,Pr/100e3,XA)] # mol/s

    # first block is region A -> just reactions
    RHS = (DE_Model.spec_MW[1:end-1]*1e-3).*[-r[1]+r[2];# production of CO in kg/s
            -r[2]-r[3];# production of CO2
            -2*r[1]-r[2]-3*r[3];# production of H2
            r[1]+r[3];# production of CH3OH
            r[2]+r[3]] # production of H2O,
    du[3:end] .= M\(RHS./(rho_STP*GHSV_end))
    du[1] = gas.density # this is the equation for getting residence time
    du[2] = 0.0
end

fluid_props.setTPX(gas,(u0[2], u0[1], u0[3:end]), spec_ind)
u0_PFR_a=[0;gas.h;u0[3:end]]

prob4=ODE.ODEProblem(PF_reactor_adiabatic!,u0_PFR_a,(0,1),saveat=zstar)
sol4 = ODE.solve(prob4,ODE.RadauIIA5(autodiff=false),
    p=param_PFR, reltol=1e-8, abstol=1e-10)
y_PFR_a=permutedims(stack(sol4.u))
GHSV_PFR_a=GHSV_constV_a[end]*y_PFR_a[end,1] ./ y_PFR_a[:,1]
Xall_PFR_a=[y_PFR_a[:,3:end] (1 .- sum(y_PFR_a[:,3:end], dims=2))]

T_PFR_a=similar(GHSV_PFR_a)
for i in eachindex(GHSV_PFR_a)
    fluid_props.setHPX(gas, (y_PFR_a[i,2], Pr, y_PFR_a[i,3:end]), spec_ind)
    T_PFR_a[i]=gas.T
end


begin
plt_PFR=plot(GHSV_PFR_a*1e3,Xall_PFR_a, label=permutedims(specs), 
    color=permutedims(cols[1:length(specs)]))
plot!(plt_PFR,GHSV_PFR_i*1e3,Xall_PFR_i, label=:none, 
    color=permutedims(cols[1:length(specs)]),
    linestyle=:dot)
plot!(plt_PFR,xaxis=:log, xlim=(1e3,1e6),ylabel="Mole Fraction",xlabel="GHSV (scc/gcat/hr)",
    title="Plug Flow Reactor\n Solid=Adiabatic, Dash=Isothermal", ylim=(0,0.30))
plot!(plt_PFR,[NaN],[NaN], label="T (ad)", color=:black, linestyle=:solid)
plot!(twinx(plt_PFR), GHSV_PFR_a*1e3, T_PFR_a .- 273.15, ylabel="T (degC)",
    color=:black, linestyle=:solid, xlim=(1e3,1e6), legend=:none, xaxis=:log, ylim=(240, 300))
end

plot(plt_constV,plt_PFR, layout=grid(1,2), size=(1080,720), margin=10mm)


##############################################
# Now we do well-stirred reactor

function WSR_reactor_adiabatic!(du, u, p)
    gas=p[1]
    uloc=copy(u)
    #uloc[u.<=0.0] .*= -1.0
    #uloc[u.>=1.0] .= 1.0
    fluid_props.setTPX(gas, (u[1], Pr, uloc[2:end-1]), spec_ind)
    T=gas.T
    h=gas.enthalpy_mass
    Y=gas.Y[spec_ind[1:end-1]]
    Yin=p[2]
    rho_STP=p[3]
    GHSV=p[4]
    hin=p[5]
    XA=uloc[2:end]
    XA[XA.<=0.0] .= 0.0
    # assume that catalyst is at wall temperature
    r=[0.0, 0.0, 0.0]

    r.= [fluid_props.r1(T,Pr/100e3,XA),
        fluid_props.r2(T,Pr/100e3,XA), 
        fluid_props.r3(T,Pr/100e3,XA)] # mol/s

    # first block is region A -> just reactions
    RHS = (DE_Model.spec_MW[1:end-1]*1e-3).*[-r[1]+r[2];# production of CO in kg/s
            -r[2]-r[3];# production of CO2
            -2*r[1]-r[2]-3*r[3];# production of H2
            r[1]+r[3];# production of CH3OH
            r[2]+r[3]] # production of H2O,
    du[2:end-1] .= (Y .- RHS./(rho_STP*GHSV) .- Yin)
    du[end] = 1.0 - sum(uloc[2:end])
    du[1]=(h-hin)*1e-3
end

u0_WSR_a=[Tr; Xall_constV_a[50,:]]
param_WSR=(gas, Yin, ρSTP, GHSV_WSR[50]/3600, hin)
fun_a(du,u)=WSR_reactor_adiabatic!(du,u,param_WSR)
sol_a=NLsolve.nlsolve(fun_a, u0_WSR_a, iterations=1000)



function WSR_reactor_isothermal!(du, u, p)
    gas=p[1]
    uloc=copy(u)
    #uloc[u.<=0.0] .*= -1.0
    #uloc[u.>=1.0] .= 1.0
    fluid_props.setTPX(gas, (Tr, Pr, uloc[1:end-1]), spec_ind)
    T=gas.T
    Y=gas.Y[spec_ind[1:end-1]]
    Yin=p[2]
    rho_STP=p[3]
    GHSV=p[4]
    XA=uloc
    XA[XA.<=0.0] .= 0.0
    # assume that catalyst is at wall temperature
    r=[0.0, 0.0, 0.0]

    r.= [fluid_props.r1(T,Pr/100e3,XA),
        fluid_props.r2(T,Pr/100e3,XA), 
        fluid_props.r3(T,Pr/100e3,XA)] # mol/s

    # first block is region A -> just reactions
    RHS = (DE_Model.spec_MW[1:end-1]*1e-3).*[-r[1]+r[2];# production of CO in kg/s
            -r[2]-r[3];# production of CO2
            -2*r[1]-r[2]-3*r[3];# production of H2
            r[1]+r[3];# production of CH3OH
            r[2]+r[3]] # production of H2O,
    du[1:end-1] .= (Y .- RHS./(rho_STP*GHSV) .- Yin)
    du[end] = 1.0 - sum(uloc)
end

using NLsolve
fluid_props.setTPX(gas,(u0[2], u0[1], u0[3:end]), spec_ind)
Yin=copy(gas.Y[spec_ind[1:end-1]])
hin=copy(gas.enthalpy_mass)
GHSV_WSR=GHSV_constV_i[2:5:end]
Xall_WSR_i=Matrix{Float64}(undef, length(GHSV_WSR), length(specs))
Xall_WSR_a=similar(Xall_WSR_i)
T_WSR_a=similar(GHSV_WSR)

for i in axes(Xall_WSR_i,1)
    @printf("Solving WSR i=%u for GHSV = %.2f\n",i, GHSV_WSR[i])
    param_WSR=(gas, Yin, ρSTP, GHSV_WSR[i]/3600, hin)
    #u0_WSR=copy(gas.X[spec_ind])
    if i>1
        u0_WSR_i=Xall_WSR_i[i-1,:]
        u0_WSR_a=[T_WSR_a[i-1]+10; Xall_WSR_a[i-1,:]]
    else
        u0_WSR_i=Xall_constV_i[2,:]
        u0_WSR_a=[Tr; Xall_constV_a[2,:]]
    end
    fun_i(du,u)=WSR_reactor_isothermal!(du,u,param_WSR)
    fun_a(du,u)=WSR_reactor_adiabatic!(du,u,param_WSR)
    sol_i=NLsolve.nlsolve(fun_i, u0_WSR_i, iterations=1000)
    Xall_WSR_i[i,:]=sol_i.zero
    sol_a=NLsolve.nlsolve(fun_a, u0_WSR_a, iterations=1000)
    Xall_WSR_a[i,:]=sol_a.zero[2:end]
    T_WSR_a[i]=sol_a.zero[1]
end



begin
plt_WSR=plot(GHSV_WSR*1e3,Xall_WSR_a, label=permutedims(specs), 
    color=permutedims(cols[1:length(specs)]))
plot!(plt_WSR,GHSV_WSR*1e3,Xall_WSR_i, label=:none, 
    color=permutedims(cols[1:length(specs)]),
    linestyle=:dot)
plot!(plt_WSR,xaxis=:log, xlim=(1e3,1e6),ylabel="Mole Fraction",xlabel="GHSV (scc/gcat/hr)",
    title="Well-Stirred Reactor\n Solid=Adiabatic, Dash=Isothermal", ylim=(0,0.30))
plot!(plt_WSR,[NaN],[NaN], label="T (ad)", color=:black, linestyle=:solid)
plot!(twinx(plt_WSR), GHSV_WSR*1e3, T_WSR_a .- 273.15, ylabel="T (degC)",
    color=:black, linestyle=:solid, xlim=(1e3,1e6), legend=:none, xaxis=:log, ylim=(240, 300))
#plot!(plt_WSR, margin=10mm)
end

###################################################################3
# plot methanol on all plots
plt_comm=Plots.plot()
plot!(plt_comm, GHSV_constV_i*1e3, Xall_constV_i[:,5], label="Const V Isothermal")
plot!(plt_comm, GHSV_constV_a*1e3, Xall_constV_a[:,5], label="Const V Adiabatic")
plot!(plt_comm, GHSV_PFR_i*1e3, Xall_PFR_i[:,5], label="PFR Isothermal")
plot!(plt_comm, GHSV_PFR_a*1e3, Xall_PFR_a[:,5], label="PFR Adiabatic")
plot!(plt_comm, GHSV_WSR*1e3, Xall_WSR_i[:,5], label="WSR Isothermal")
plot!(plt_comm, GHSV_WSR*1e3, Xall_WSR_a[:,5], label="WSR Adiabatic")
plot!(plt_comm, xaxis=:log, xlim=(1e3,1e6),ylabel="Mole Fraction CH3OH",xlabel="GHSV (scc/gcat/hr)",
    title="Comparison of Methanol Production", ylim=(0,0.15), legend=:topright)






fig=plot(plt_constV,plt_PFR,plt_WSR, plt_comm, layout=grid(2,2),
size=(1440,1080), margin=10mm)
Tstr=@sprintf("%.0f",Tr - 273.15)
savefig(fig, topdir*"Misc/Representative_Figures/Reactor_Comparison_"*Tstr*"C.png")
