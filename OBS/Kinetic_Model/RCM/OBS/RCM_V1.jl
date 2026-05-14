topdir=abspath((@__DIR__)*"/../")
using Pkg
Pkg.activate(topdir*"/src/.")

using LinearAlgebra

include(topdir*"src/fluid_properties/fluid_props.jl")


gas1=fluid_props.initialize_ideal_gas("gri30.yaml")
gas2=fluid_props.initialize_ideal_gas("gri30.yaml")
fluid_props.kingas=fluid_props.initialize_ideal_gas("gri30.yaml")
rhocat=1e3 # kg/m^3


specs=["CO", "CO2", "H2", "CH3OH", "H2O", "N2"]

# index of species in cantera object (+1 to convert to julia indexing)
spec_ind=fluid_props.spec_inds(gas1,specs) .+ 1
spec_MW=fluid_props.MW_spec(gas1)[spec_ind]
Nspec=length(specs)
u0=fluid_props.u0(gas1, specs) # ref. state internal energy of each species (J/kg)
u0_mol=u0 .* spec_MW /1000 # ref. state internal energy (J/mol)


function Mmat(gasA, gasB, Vstar)
    return [fluid_props.int_nrg(gasA)-fluid_props.enthalpy(gasB) 0 fluid_props.rho(gasA) 0;
        fluid_props.enthalpy(gasB)  Vstar*fluid_props.int_nrg(gasB) 0 Vstar*fluid_props.rho(gasB);
        1   Vstar   0   0]
end

function Dmat(gasA, gasB)
    return [fluid_props.rho(gasA)/fluid_props.P(gasA) -fluid_props.rho(gasA)/fluid_props.T(gasA)    0;
        fluid_props.rho(gasB)/fluid_props.P(gasB)   0   -fluid_props.rho(gasB)/fluid_props.T(gasB);
        0   fluid_props.cv(gasA)    0;
        0   0   fluid_props.cv(gasB)]
end

function Bmat(gas)
    MW=fluid_props.MW_mix(gas)
    return 1/MW * (LinearAlgebra.diagm(spec_MW[1:end-1]) - 
        1/MW*(transpose(spec_MW[1:end-1] .- spec_MW[end])) .* fill(1.0,(length(spec_MW)-1,1))) 
end



function Mfull(gasA, gasB, V)
    # M is the coefficient matrix of the base equation
    M=zeros(3+2*(Nspec-1),4+2*(Nspec-1))
    M[1:3,1:4] .= Mmat(gasA, gasB, V)
    YA=fluid_props.Y(gasA)[spec_ind]
    YB=fluid_props.Y(gasB)[spec_ind]
    rhoA=fluid_props.rho(gasA)
    rhoB=fluid_props.rho(gasB)
    for i=1:Nspec-1
        M[3+i,1]=YA[i]-YB[i]
        M[3+i, 4+i] = rhoA
        M[3+Nspec-1+i, 1] = YB[i]
        M[3+Nspec-1+i, 2] = YB[i]*V
        M[3+Nspec-1+i, 4+Nspec-1+i] = rhoB*V
    end
    return M
end

function Dfull(gasA, gasB)
    # D converts rho and u to P and T
    D=zeros(4+2*(Nspec-1), 3+2*(Nspec-1))
    D[1:4, 1:3] .= Dmat(gasA, gasB)
    D[5:end,4:end] .= LinearAlgebra.diagm(ones(2*(Nspec-1)))
    for i=1:Nspec-1
        D[1,3+i] = 1/fluid_props.MW_mix(gasA)*(spec_MW[i]-spec_MW[end])
        D[2,3+Nspec-1+i] = 1/fluid_props.MW_mix(gasB)*(spec_MW[i]-spec_MW[end])
    end
    return D
end

function Bfull(gasA, gasB)
    # B converts mass fraction to mole fraction
    B=LinearAlgebra.diagm(ones(4+2*(Nspec-1)))
    B[5:5+(Nspec-2),5:5+(Nspec-2)] .= Bmat(gasA)
    B[5+Nspec-1:end,5+Nspec-1:end] .= Bmat(gasB)
    return B
end


fluid_props.setTPX(gas1, (300.0, 100e3, "H2:0.6, CO:0.25, CO2:0.05, N2:0.1"))
fluid_props.setTPX(gas2, (300.0, 100e3, "H2:0.6, CO:0.25, CO2:0.05, N2:0.1"))


M=Mfull(gas1, gas2, 10)
D=Dfull(gas1,gas2)
B=Bfull(gas1,gas2)

function V_dV(t, CompRatio, tcomp)
    dV=-(CompRatio-1)/tcomp
    if t<tcomp-1e-5
        return (CompRatio-1) + dV*t, dV
    else
        return (CompRatio-1) + dV*(tcomp-1e-5), 0.0
    end 
end

# Non-reactive
function f!(du, u, p, t)
    gasA=p[1] # for passing to fluid props
    gasB=p[2] # for passing to fluid props
    V, dV = p[3](t) # function of t returning V and dV
    fluid_props.setTPX(gasA, (u[2], u[1], gasA.X))
    fluid_props.setTPX(gasB, (u[3], u[1], gasB.X))
    M=Mfull(gasA, gasB, V)
    D=Dfull(gasA, gasB)
    B=Bfull(gasA, gasB)
    M2=(M*B*D)
    du .= M2\vcat([0; -fluid_props.rho(gasB)*fluid_props.enthalpy(gasB)*dV; -fluid_props.rho(gasB)*dV],
        zeros(Nspec-1),
        -fluid_props.Y(gasB)[spec_ind[1:end-1]].*fluid_props.rho(gasB)*dV)
    return nothing
end


CR=16
Vfunc(t)=V_dV(t, CR, 10e-3)
fluid_props.setTPX(gas1, (300.0, 100e3, "H2:0.6, CO:0.25, CO2:0.05, N2:0.1"))
fluid_props.setTPX(gas2, (300.0, 100e3, "H2:0.6, CO:0.25, CO2:0.05, N2:0.1"))
param=(gas1, gas2, Vfunc)

u0=vcat([100e3, 300.0, 300.0], 
    fluid_props.X(gas1)[spec_ind[1:end-1]], 
    fluid_props.X(gas2)[spec_ind[1:end-1]])
du=copy(u0)

f!(du, u0, param, 0.0)
using OrdinaryDiffEq
tmax=8e-3;
prob=ODEProblem(f!,u0,(0.0,tmax))
sol = solve(prob,Tsit5(),p=param)

tp=range(0.0,tmax,101)
y=permutedims(stack(sol.(tp)))
V=permutedims(stack(Vfunc.(tp)))[:,1]

fluid_props.setTPX(gas1, (300.0,100e3,gas1.X))
gamma=fluid_props.cp(gas1)/fluid_props.cv(gas1)

using Plots
using LaTeXStrings
begin
fig=plot(tp*1e3,y[:,1]/100e3, label=L"\frac{P_2}{P_1}", xlabel="Time (ms)")
plot!(tp*1e3,( (1 .+ V[1])./(1 .+ V)).^gamma, label=L"\left(\frac{V_1}{V_2}\right)^\gamma")
plot!(tp*1e3,( y[:,2]./300.0).^(gamma/(gamma-1)), label=L"\left(\frac{T_2^A}{T_1^A}\right)^\frac{\gamma}{\gamma-1}")
plot!(tp*1e3,( y[:,3]./300.0).^(gamma/(gamma-1)), label=L"\left(\frac{T_2^B}{T_1^B}\right)^\frac{\gamma}{\gamma-1}")#yscale=:log)
end
savefig(fig,topdir*"output/figures/Isentropic_Compression_demo.png")

x=(1 .+ V[1])./(1 .+ V)
begin
    fig2=plot(x, 100e3*x.^gamma, xscale=:log, yscale=:log)
    plot!(x, y[:,1])
    plot!(x,100e3*( y[:,2]./300.0).^(gamma/(gamma-1)))
    plot!(x,100e3*( y[:,3]./300.0).^(gamma/(gamma-1)))
end
savefig(fig2,topdir*"output/figures/Isentropic_Compression_demo2.png")


###################################################
# Reactive

function V_dV(t, CompRatio, tcomp)
    dV=-(CompRatio-1)/tcomp
    if t<tcomp
        return (CompRatio-1) + dV*t, dV
    else
        return 0.0, 0.0
    end 
end

function f2!(du, u, p, t)
    gasA=p[1] # for passing to fluid props
    gasB=p[2] # for passing to fluid props
    V, dV = p[3](t) # function of t returning V and dV
    XA=0.0*gasA.X
    XA[spec_ind[1:end-1]].=u[4:4+length(spec_ind)-2]
    XA[spec_ind[end]] = 1-sum(XA)
    fluid_props.setTPX(gasA, (u[2], u[1], XA))
    XB=0.0*gasB.X;
    XB[spec_ind[1:end-1]].=u[4+length(spec_ind)-1:end]
    XB[spec_ind[end]] = 1-sum(XB)
    fluid_props.setTPX(gasB, (u[3], u[1], XB))
    M=Mfull(gasA, gasB, V)
    D=Dfull(gasA, gasB)
    B=Bfull(gasA, gasB)

    M2=(M*B*D)
    XA=fluid_props.X(gasA)[spec_ind]
    r= rhocat*[fluid_props.r1(u[2],u[1]/100e3,XA),
        fluid_props.r2(u[2],u[1]/100e3,XA), 
        fluid_props.r3(u[2],u[1]/100e3,XA)] # mol/s

    du .= M2\vcat([-1e5*(u[2]-(250+273.15)); -fluid_props.rho(gasB)*fluid_props.enthalpy(gasB)*dV; -fluid_props.rho(gasB)*dV],
        (spec_MW[1:end-1]*1e-3).*[-r[1]+r[2];# production of CO in kg/s
            -r[2]-r[3];# production of CO2
            -2*r[1]-r[2]-3*r[3];# production of H2
            r[1]+r[3];# production of CH3OH
            r[2]+r[3]],# production of H2O,
        -fluid_props.Y(gasB)[spec_ind[1:end-1]].*fluid_props.rho(gasB)*dV)

    return nothing
end


CR=12
tcomp=50e-3
Vfunc(t)=V_dV(t, CR, tcomp)
fluid_props.setTPX(gas1, (300.0, 100e3, "H2:0.6, CO:0.25, CO2:0.05, N2:0.1"))
fluid_props.setTPX(gas2, (300.0, 100e3, "H2:0.6, CO:0.25, CO2:0.05, N2:0.1"))
param=(gas1, gas2, Vfunc)

u0=vcat([5*100e3, 300.0, 300.0], 
    fluid_props.X(gas1)[spec_ind[1:end-1]], 
    fluid_props.X(gas2)[spec_ind[1:end-1]])
du=copy(u0)


tmax=tcomp;
prob=ODEProblem(f2!,u0,(0.0,tmax))
#sol = solve(prob,Tsit5(),p=param)
sol = solve(prob,Rosenbrock23(autodiff=false),p=param)

tp=range(0.0,tmax,101)
y=permutedims(stack(sol.(tp)))
V=permutedims(stack(Vfunc.(tp)))[:,1]

fluid_props.setTPX(gas1, (300.0,100e3,gas1.X))
gamma=fluid_props.cp(gas1)/fluid_props.cv(gas1)



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
savefig(fig,topdir*"output/figures/rxn_prelim2")