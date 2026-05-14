#################################3
# Modelled using the reaction kinetics from Portha et al. 2017
# Ind. Eng. Chem. Res. 2017, 56, 13133-13145

cd(@__DIR__)
using Pkg
Pkg.activate(".")

using PyCall
using Printf
using Plots
using OrdinaryDiffEq
ct=pyimport("cantera")
gas=ct.Solution("gri30.yaml")
const Ru=8.314 # J/mol K

function G(T,P,spec)
    gas.TPX=T,P, @sprintf("%s:1", spec)
    return gas.enthalpy_mole - gas.entropy_mole * T
end

# equilibrium constants
function K_CO2hyd(T,P)
    dG=G(T,P,"CH3OH")+G(T,P,"H2O")-
        G(T,P,"CO2")-3*G(T,P,"H2")
    return exp(-dG/(Ru*1e3*T))*(100e3/P)^2 # bar^-2
end
function K_COhyd(T,P)
    dG=G(T,P,"CH3OH")-
        G(T,P,"CO")-2*G(T,P,"H2")
    return exp(-dG/(Ru*1e3*T))*(100e3/P)^2 # bar^-2
end
function K_RWGS(T,P)
    dG=G(T,P,"CO2")+G(T,P,"H2") - 
    G(T,P,"CO")-G(T,P,"H2O")
    return exp(-dG/(Ru*1e3*T))
end

# adsorption constants (Table 7)
bCO(T)=2.16e-5*exp(46800/(Ru*T)) # bar ^-1
bCO2(T)=7.05e-7*exp(61700/(Ru*T)) # bar ^-1
bH2O_H2(T)=6.37e-9*exp(84000/(Ru*T)) # bar ^-1/2

# kinetic rates (Table 8)
k_COhyd(T) = 4.89e7*exp(-113000/(Ru*T)) # mol/s/bar/kgcat
k_RWGS(T) = 9.64e11*exp(-152900/(Ru*T)) # mol/s/bar^(1/2)/kgcat
k_CO2hyd(T) = 1.09e5*exp(-87500/(Ru*T)) # mol/s/bar/kgcat

# reaction rates
# x will be mole fractions ordered as [CO, CO2, H2, CH3OH, H2O]
denom(T, P, x)=(1 .+ bCO(T)*P*x[1]+bCO2(T)*P*x[2]).*
    ((P*x[3]).^(1/2)+(bH2O_H2(T))*P*x[5])

r1(T,P,x) = k_COhyd(T)*bCO(T)*(
    (P*x[1].*(P*x[3]).^(3/2)-P*x[4]./((P*x[3]).^(1/2)*K_COhyd(T,P)))./
    denom(T,P,x))
r2(T,P,x) = k_RWGS(T)*bCO2(T)*(
        (P*x[2].*P.*x[3]-(P*x[1].*P.*x[5])./(K_RWGS(T,P)))./
        denom(T,P,x))
r3(T,P,x) = k_CO2hyd(T)*bCO2(T)*(
    (P*x[2].*(P*x[3]).^(3/2)-P*x[4]*P.*x[5]./((P*x[3]).^(3/2)*K_CO2hyd(T,P)))./
    denom(T,P,x))

rho_c=135e-6/(0.307e-6) # density of catalyst

function f!(du,u,p,t)
    # u is the vector of molar densities
    # p (parameters) is [Ptot, Ttot]
    T=p[1];
    nt=sum(u)
    P=nt*Ru*T
    x=u/nt
    r=rho_c*[r1(T,P,x),r2(T,P,x), r3(T,P,x)]
    du[1]=-r[1]+r[2]# production of CO
    du[2]=-r[2]-r[3]# production of CO2
    du[3]=-2*r[1]-r[2]-3*r[3]# production of H2
    du[4]=r[1]+r[3]# production of CH3OH
    du[5]=r[2]+r[3]# production of H2O
end

P0=200e5
T0=500.0
x0=[0.0, 0.25, 0.75, 0.0,0.0]
u0=x0/sum(x0) * P0/(Ru*T0)
du0=zeros(Float64, 5)
prob=ODEProblem(f!, u0, (0.0,0.5))
sol=solve(prob,Tsit5(), p=(T0))

# simple version - No CO reactions
function fa!(du,u,p,t)
    # u is the vector of molar densities
    # p (parameters) is [Ptot, Ttot]
    pp=u * Ru * p[1] # partial pressures of all species
    rf=4.9e-7*exp(-64e3/(Ru*p[1]))*
        pp[1]^(0.23)*pp[2]^1.7 *
        rho_c# rate of CH3OH production (unsure units)
    du[1]=-rf# production of CO2
    du[2]=-3*rf# production of H2
    du[3]=rf# production of CH3OH
    du[4]=rf# production of H2
end

u0a=[0.25, 0.75, 0,0] * P0/(Ru*T0)
du0a=zeros(Float64, 4)
proba=ODEProblem(fa!, u0a, (0.0,0.5))
sola=solve(proba,Tsit5(),p=(T0))



lab=["CO","CO2","H2","CH3OH","H2O"]
tsamp=range(0,1,501)
n=permutedims(hcat(map(t -> sol(t), tsamp)...))
nt=sum(n, dims=2)
na=permutedims(hcat(map(t -> sola(t), tsamp)...))
nta=sum(na, dims=2)
Conv=(u0[2] .- n[:,2])/u0[2]
Conva=(u0a[1] .- na[:,1])/u0a[1]

begin
    plot(tsamp*1e3, n ./ nt, label=permutedims(lab), 
        xlabel="Time (ms)", ylabel="Mole Fraction")
    plot!(tsamp*1e3, na ./ nta, label=permutedims(lab[2:end]), 
        xlabel="Time (ms)", ylabel="Mole Fraction", linestyle=:dash)
    ax_right=plot!(twinx(), tsamp*1e3, Conv, color=:black, 
        label="conversion", ylim=(0,1))
    plot!(twinx(),tsamp*1e3, Conva, color=:black, 
        linestyle=:dash, ylim=(0,1));
end


rho_in=[28, 44, 2, 32, 18] .* u0



begin
ax=Plots.plot()
Prng_bar=vcat(collect(4:4:36),collect(40:40:400))
for P0=Prng_bar*1e5
T0=580.0
x0=[0.0, 0.25, 0.75, 0.0,0.0]
u0=x0/sum(x0) * P0/(Ru*T0)
du0=zeros(Float64, 5)
prob=ODEProblem(f!, u0, (0.0,1.0))
tsamp=range(0,1,1001)
sol=solve(prob,Tsit5(), p=(T0))
n=permutedims(hcat(map(t -> sol(t), tsamp)...))
nt=sum(n, dims=2)
Conv=(u0[2] .- n[:,2])/u0[2]
Vin_actual=sum(u0)*Ru*T0/P0 # assume this is the catalyst volume
Vin_STP=sum(u0)*Ru*298.15/100e3 # STP volume of gas into reactor
GHSV=Vin_STP/Vin_actual * 1 ./(tsamp) * 3600
lab=@sprintf("P=%.0f bar", P0/100e3)
plot!(ax,GHSV, Conv*100, xaxis=:log,
    xlabel="GHSV (hr^-1)", ylabel="CO2 Conversion (%)", label=lab)
end
plot(ax)
end