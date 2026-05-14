topdir=abspath((@__DIR__)*"./../../")
using Pkg
Pkg.activate(topdir*"Kinetic_Model/.")

import OrdinaryDiffEq  as ODE

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
Treactor=290+273.15
comp="CO2:0.25, H2:0.75";

# Get STP conditions
gas.TPX=(273.15, 101.325e3, comp)
ρSTP=gas.density

GHSV_min=8 #L/gcat/hr
#################################################
# Now do PFR
# use formulation I alread developed (python model)
zstar=range(0, 1, 201)

# isothermal
function PF_reactor_isothermal!(du, u, p, t)
    gas=p[1]
    Tr=p[5]
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

# function for solving PFR
function solve_PFR(comp, GHSV, Tr)
    gas.TPX=(Tr, Pr, comp)
    u0=[0; fluid_props.X(gas)[spec_ind[1:end-1]]]
    param_PFR=(gas, rhocat_bulk, ρSTP, GHSV/3600, Tr)
    prob=ODE.ODEProblem(PF_reactor_isothermal!,u0,(0,1),saveat=zstar)
    sol = ODE.solve(prob,ODE.RadauIIA5(autodiff=false),
        p=param_PFR, reltol=1e-8, abstol=1e-10)
    y=permutedims(stack(sol.u))
    G_inst=GHSV_min*y[end,1] ./ y[:,1]
    X=[y[:,2:end] (1 .- sum(y[:,2:end], dims=2))]
    return G_inst, X
end

GHSV, X=solve_PFR(comp, 8, 273.15+250.0)

using Plots, Measures, Printf, LaTeXStrings, ColorSchemes
cols=ColorSchemes.tab10
begin
plt=plot(GHSV*1e3,X, label=permutedims(specs), 
    color=permutedims(cols[1:length(specs)]))
plot!(plt,xaxis=:log, xlim=(6e3,1e6),ylabel="Mole Fraction",xlabel="GHSV (scc/gcat/hr)",
    ylim=(0,0.30))
end

##############################################

vec2compstr(vec)=join([string(specs[i], ":", vec[i]) for i in eachindex(specs)], ", ")
MW_mix(X)=sum(X.*DE_Model.spec_MW)

function solve_recirc!(du,u,SN_target, GHSV, Tr, purge_frac=0.01)
    # u is the reactor inlet composition (mole fractions)
    # first solve the PFR
    comp_str=vec2compstr(u)
    MW_in=MW_mix(u)
    _,X=solve_PFR(comp_str, GHSV, Tr)
    # for now assume perfect separation
    Xout=X[end,:]
    MW_out=MW_mix(Xout)
    # moles out (per mole in)
    nout=MW_in/MW_out
    Xout[[4,5]] .= 0.0
    nrecirc=nout*(sum(Xout))
    nsep=nout-nrecirc
    nrecirc *= (1-purge_frac)
    Xout .*= nrecirc/sum(Xout)
    # now we need to correct SN
    SN0=(Xout[3]-Xout[2])/(Xout[1]+Xout[2])
    if SN0<SN_target
        # we need to add H2
        Xout[3]=(Xout[1]+Xout[2])*SN_target + Xout[2]
    else
        # we need to add CO2 (unlikely)
        Xout[2]=-(SN_target*Xout[1]-Xout[3])/(SN_target+1)
    end
    # finally, we add enough CO2 and H2 (at the target SN) to get 1 mole total
    nadd=1-sum(Xout)
    Xout[2] += nadd/(2+SN_target)
    Xout[3] += nadd-nadd/(2+SN_target)
    du .= Xout - u
end


u0=X[end,:]
u0[[4,5]] .= 0.0
u0 ./= sum(u0)
using NLsolve
du=copy(u0)
SN_target=2.0
fun!(du,u) = solve_recirc!(du,u,SN_target, 20, 273.15+270)
res=NLsolve.nlsolve(fun!, u0, iterations=60,show_trace=true)

fun!(du,u0)


G_list=10 .^ range(1,3,21)
T_list=273.15 .+ range(250,300,6)

X=Array{Float64}(undef, length(G_list), length(specs), length(T_list))
Xin=similar(X)
for (i, G) in enumerate(G_list), (j, T) in enumerate(T_list)
    @printf("Solving for GHSV = %d, T=%d...", G, T)
    fun!(du,u) = solve_recirc!(du,u,SN_target, G, T, 0.02)
    res=NLsolve.nlsolve(fun!, u0, iterations=40)
    _, Xtmp = solve_PFR(vec2compstr(res.zero), G,T)
    X[i,:,j] = Xtmp[end,:]
    Xin[i,:,j] = res.zero
    u0=res.zero
    res.f_converged ? print("Success\n") : print("Failed\n")
end



begin
    plt=[Plots.plot() for _ in 1:6]
    for (i,p) in enumerate(plt)
        contourf!(p, G_list*1e3, T_list .-273.15, Matrix(X[:,i,:])'*100, levels=50, linewidth=0,
            title=@sprintf("Reactor Outlet %s%%", specs[i]), xaxis=:log, xlabel="GHSV (scc/gcat/hr)", ylabel="Temperature (°C)")
    end
    plot(plt..., layout=(3,2), size=(1440,1440))
end



###
#=
plot(G_list, X, label=permutedims(specs), 
    color=permutedims(cols[1:length(specs)]),
    ylim=(0,0.3))


# single pass
G, X_SP=solve_PFR(comp, G_list[1])

begin
plt=plot(G*1e3,X_SP, label=permutedims(specs), 
    color=permutedims(cols[1:length(specs)]))
plot!(plt, G_list*1e3, X, label=:none, color=permutedims(cols[1:length(specs)]), linestyle=:dash)
plot!(plt,xaxis=:log, xlim=(10e3,1e6),ylabel="Mole Fraction",xlabel="GHSV (scc/gcat/hr)",
    ylim=(0,0.30), title="Reactor Outlet\n Solid=SP, Dashed=Recirc")
end
begin
plt=plot(G_list*1e3,Xin, label=permutedims(specs), 
    color=permutedims(cols[1:length(specs)]))
plot!(plt, G_list*1e3, X, label=:none, color=permutedims(cols[1:length(specs)]), linestyle=:dash)
plot!(plt,xaxis=:log, xlim=(10e3,1e6),ylabel="Mole Fraction",xlabel="GHSV (scc/gcat/hr)",
    ylim=(0,0.80), title="Reactor Outlet\n Solid=Inlet, Dashed=Outlet")
end
=#