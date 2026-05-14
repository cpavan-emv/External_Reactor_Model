#################################3
# Modelled using the reaction kinetics from Portha et al. 2017
# Ind. Eng. Chem. Res. 2017, 56, 13133-13145

topdir=abspath((@__DIR__)*"./../")
using Pkg
Pkg.activate(topdir*"/src/.")


using PyCall
using Printf
using Plots
using OrdinaryDiffEq
using LinearAlgebra
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
# units should be mol/s/kgcat when P is in bar
denom(T, P, x)=(1 .+ bCO(T)*P*x[1]+bCO2(T)*P*x[2]).*
    ((P*x[3]).^(1/2)+(bH2O_H2(T))*P*x[5])

r1(T,P,x) = k_COhyd(T)*bCO(T)*(
    (P*x[1].*(P*x[3]).^(3/2)-P*x[4]./((P*x[3]).^(1/2)*K_COhyd(T,100e3)))./
    denom(T,P,x))
r2(T,P,x) = k_RWGS(T)*bCO2(T)*(
        (P*x[2].*P.*x[3]-(P*x[1].*P.*x[5])./(K_RWGS(T,100e3)))./
        denom(T,P,x))
r3(T,P,x) = k_CO2hyd(T)*bCO2(T)*(
    (P*x[2].*(P*x[3]).^(3/2)-P*x[4]*P.*x[5]./((P*x[3]).^(3/2)*K_CO2hyd(T,100e3)))./
    denom(T,P,x))


rho_c=3e-3/(4.5e-3^2*pi*4.5e-2) # density of catalyst (kg/m3)
MW=[28, 44, 2, 32, 18, 28] # g/mol

function MWbar(x)
    return transpose(MW) * x
end

function M(x)
    M=zeros(length(x), length(x))
    for i in eachindex(x)[1:end-1]
        M[i,:] .= -x[i] .* (MW .- MW[end]) / MWbar(x)
        M[i,i] += 1
    end
    M[:, end] .=0
    M[end,:] .= 0
    return M
end

function f!(du,u,p,t)
    # u is the vector of mole fractions
    # p (parameters) is [T, Ptot, rho0*v0]
    T=p[1]; P=p[2]
    # RHS equations
    x=u
    r=rho_c*[r1(T,P/100e3,x),r2(T,P/100e3,x), r3(T,P/100e3,x)]*MWbar(x)/p[3]   # units g/(m3*s) * 1/(g/m2*s) = 1/m
    fun=Vector{Float64}(undef, 6)
    fun[1]=-r[1]+r[2]# production of CO
    fun[2]=-r[2]-r[3]# production of CO2
    fun[3]=-2*r[1]-r[2]-3*r[3]# production of H2
    fun[4]=r[1]+r[3]# production of CH3OH
    fun[5]=r[2]+r[3]# production of H2O
    fun[6]=1-sum(x) # algebraic constraint equation gives N2 mole fraction

    du[1:5] = M(x)[1:5,1:5] \ fun[1:5]
    du[6] = 1-sum(u)
end

Mconst=diagm([ones(5);0])

Axsect=9e-3^2*pi/4
P0=50*100e3
T0=240+273.15
x0=[0.3, 0.025, 0.54, 0.0,0.0,0.135]
x0 /= sum(x0)
rho0=P0/(Ru*T0)*MWbar(x0) # g/m3
v0=100/(Axsect*1e4)/60/100 # flow speed in m/s for 100sccm flow through 0.5cm diameter channel
flux0=rho0*v0 # g/(m2*s) - > the factor of 100 needs to be here to match order of magnitude of exp, but can't find why...

du0=zeros(Float64, 6)
fun_test=ODEFunction(f!, mass_matrix=Mconst)
prob=ODEProblem(fun_test, x0, (0.0,10.0))
sol=solve(prob, Rodas4(autodiff=false),p=(T0,P0, flux0))

#prob=ODEProblem(f!, u0, (0.0,1), mass_matrix=Mconst)
#sol=solve(prob, Rodas4(), p=(T0,P0, flux0))


lab=["CO","CO2","H2","CH3OH","H2O","N2"]
zsamp=range(0,10,501)
x=permutedims(hcat(map(z -> sol(z), zsamp)...))
MW_gas=stack([MWbar(y) for y in eachrow(x)])
n=permutedims(stack([x[i,:]*MW_gas[1]/MW_gas[i] for i in eachindex(MW_gas)]))
Conv=(x0[1] .- n[:,1])/x0[1]

begin
    plot(zsamp, x, label=permutedims(lab), 
        xlabel="Position (m)", ylabel="Mole Fraction", xlim=(0,10))
    ax_right=plot!(twinx(), zsamp, Conv, color=:black, 
        label="conversion", ylim=(0,1), legend=:topright, xlim=(0,10))
end


begin
    plot(zsamp, n, label=permutedims(lab), 
        xlabel="Position (m)", ylabel="Mole Flowrate")
end


########################################################
# try without transformation into linear reactor space
function f2!(du,u,p,t)
    # u is the vector of mole fractions
    # p (parameters) is [T, Ptot, rho0*v0]
    T=p[1]; P=p[2]
    # RHS equations
    x=u / sum(u)
    r=[r1(T,P/100e3,x),r2(T,P/100e3,x), r3(T,P/100e3,x)]   # units mol/s/kgcat
    du[1]=-r[1]+r[2]# production of CO
    du[2]=-r[2]-r[3]# production of CO2
    du[3]=-2*r[1]-r[2]-3*r[3]# production of H2
    du[4]=r[1]+r[3]# production of CH3OH
    du[5]=r[2]+r[3]# production of H2O
    du[6]=0.0 # production of N2
end

ndot0 = 100/60*1e-6 * P0/(Ru*T0)
u0 = x0*ndot0

du0=zeros(Float64, 6)
fun_test2=ODEFunction(f2!)
prob2=ODEProblem(fun_test2, u0, (0.0,1))
sol2=solve(prob2, Rodas4(autodiff=false),p=(T0,P0))

msamp=range(0,1,501)
n2=permutedims(hcat(map(z -> sol2(z), msamp)...))
MW_gas=stack([MWbar(y) for y in eachrow(x)])
x2 = n2 ./ sum(n2,dims=2)
Conv2=(u0[1] .- n2[:,1])/u0[1]



L=msamp ./ rho_c /Axsect

begin
    plot(L, x2, label=permutedims(lab), 
        xlabel="Distance (m)", ylabel="Mole Fraction", xlim=(0,10))
    ax_right=plot!(twinx(), L, Conv2, color=:black, 
        label="conversion", ylim=(0,1), legend=:topright, xlim=(0,10))
end
