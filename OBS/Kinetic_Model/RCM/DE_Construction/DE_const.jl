module DE_construct
# Functions for constructing the differential equation
using LinearAlgebra
using NLsolve

include((@__DIR__)*"/DE_fully_implicit.jl")

spec_ind=[0] # index of all species
spec_MW=[0.0] # molecular weight of all species
Nspec=0 # number of species
intE0=[0.0] # ref. state internal energy of each species (J/kg)
fluid_props=nothing


function set_gas_constants(gas, specs)
    # this is basically the initializer
    DE_construct.spec_ind=fluid_props.spec_inds(gas,specs) .+ 1
    DE_construct.spec_MW=fluid_props.MW_spec(gas)[spec_ind]
    DE_construct.Nspec=length(specs)
    # can use either h or u here since everything becomes deltas
    # stick with h because I like enthalpy of formation
    #DE_construct.intE0=fluid_props.u0(gas, specs) # ref. state internal energy of each species (J/kg)
    DE_construct.intE0=fluid_props.h0(gas, specs) # ref. state internal energy of each species (J/kg)
end


function Mmat(gasA, gasB, V_dV)
    # M is the coefficients in the governing equations (excl. species cons.)
    if V_dV[2]<=0
        flow=fluid_props.enthalpy(gasB)
    elseif V_dV[2]>0
        flow=fluid_props.enthalpy(gasA)
    else
        flow=0
    end
    return [fluid_props.int_nrg(gasA)-flow 0 fluid_props.rho(gasA) 0;
        flow  V_dV[1]*fluid_props.int_nrg(gasB) 0 V_dV[1]*fluid_props.rho(gasB);
        1   V_dV[1]   0   0]
end

function Dmat(gasA, gasB)
    # D is the coefficients to convert rho to a function of P and T
    return [fluid_props.rho(gasA)/fluid_props.P(gasA) -fluid_props.rho(gasA)/fluid_props.T(gasA)    0;
        fluid_props.rho(gasB)/fluid_props.P(gasB)   0   -fluid_props.rho(gasB)/fluid_props.T(gasB);
        0   1.0    0;
        0   0   1.0]
end

function Bmat(gas)
    # B is the coefficients to convert from mass fraction to mole fraction
    MW=fluid_props.MW_mix(gas)
    X=fluid_props.X(gas)[spec_ind]
    return 1/MW * (LinearAlgebra.diagm(spec_MW[1:end-1]) - 
        1/MW*(transpose(spec_MW[1:end-1] .- spec_MW[end])) .* (spec_MW[1:end-1].*X[1:end-1])) 
end

function Mfull(gasA, gasB, V_dV)
    # M is the coefficient matrix of the base equation
    M=zeros(3+2*(Nspec-1),4+2*(Nspec-1))
    M[1:3,1:4] .= Mmat(gasA, gasB, V_dV)
    YA=fluid_props.Y(gasA)[spec_ind]
    YB=fluid_props.Y(gasB)[spec_ind]
    rhoA=fluid_props.rho(gasA)
    rhoB=fluid_props.rho(gasB)

    if V_dV[2]<=0
        case=1
    elseif V_dV[2]>0
        case=2
    else
        case=0
    end

    for i=1:Nspec-1
        M[3+i, 4+i] = rhoA
        if case==1
            M[3+i,1]=YA[i]-YB[i]
            M[3+Nspec-1+i, 1] = YB[i]      
        elseif case==2
            M[3+i,1]=0.0
            M[3+Nspec-1+i, 1] = YA[i]
        end
        M[3+Nspec-1+i, 2] = YB[i]*V_dV[1]
        M[3+Nspec-1+i, 4+Nspec-1+i] = rhoB*V_dV[1]
    end
    return M
end

function Dfull(gasA, gasB)
    # D converts rho to P and T
    D=zeros(4+2*(Nspec-1), 3+2*(Nspec-1))
    D[1:4, 1:3] .= Dmat(gasA, gasB)
    D[5:end,4:end] .= LinearAlgebra.diagm(ones(2*(Nspec-1)))
    for i=1:Nspec-1
        D[1,3+i] = fluid_props.rho(gasA)/fluid_props.MW_mix(gasA)*(spec_MW[i]-spec_MW[end])
        D[2,3+Nspec-1+i] = fluid_props.rho(gasB)/fluid_props.MW_mix(gasB)*(spec_MW[i]-spec_MW[end])
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

function Ffull(gasA, gasB)
    # F converts du to cv*dt + sum(yi * u0_i)
    F=LinearAlgebra.diagm(ones(4+2*(Nspec-1)))
    F[3,3]=fluid_props.cv(gasA)
    F[4,4]=fluid_props.cv(gasB)
    for i=1:Nspec-1
        F[3,4+i] = (intE0[i]-intE0[end])
        F[4,4+Nspec-1+i] = (intE0[i]-intE0[end])
    end
    return F
end

function setMassMat(u,p,t)
    gasA=p[1] # for passing to fluid props
    gasB=p[2] # for passing to fluid props
    V_dV = p[3](t) # function of t returning V and dV
    fluid_props.setTPX(gasA, (u[2], u[1], u[4:4+length(spec_ind)-2]), spec_ind)
    fluid_props.setTPX(gasB, (u[3], u[1],  u[4+length(spec_ind)-1:end]), spec_ind)

    # inverting before storing is only worth it if the matrix is used multiple times
    # -> potential for split problem with mass matrix cached
    # caching in split problem is tough without knowing what happens under the hood
    # would need to check for any changes to state each time

    # return LinearAlgebra.lu(Mfull(gasA, gasB, V)*
    #     Ffull(gasA,gasB)*
    #     Bfull(gasA,gasB)*
    #     Dfull(gasA,gasB))
    return Mfull(gasA, gasB, V_dV)*
         Ffull(gasA,gasB)*
         Bfull(gasA,gasB)*
         Dfull(gasA,gasB)
end

function f!(du, u, p, t)
    # this is the entire RHS of the function
    # should only be used with stiff equation solvers
    Mass=setMassMat(u,p,t)
    gasA=p[1] # for passing to fluid props
    gasB=p[2] # for passing to fluid props
    V, dV = p[3](t) # function of t returning V and dV
    rhocat=p[4]
    Tcat=p[5]
    #Tcat=250

    XA=fluid_props.X(p[1])[spec_ind]
    # r= rhocat*[fluid_props.r1(u[2],u[1]/100e3,XA),
    #     fluid_props.r2(u[2],u[1]/100e3,XA), 
    #     fluid_props.r3(u[2],u[1]/100e3,XA)] # mol/s
    # # fix catalyst temp
    r= rhocat*[fluid_props.r1(Tcat,u[1]/100e3,XA),
        fluid_props.r2(Tcat,u[1]/100e3,XA), 
        fluid_props.r3(Tcat,u[1]/100e3,XA)] # mol/s

    τ=50/1000#τ_lookup(t)
    k=fluid_props.rho(gasA)*fluid_props.cp(gasA)/τ
    #println(k)


    du .= Mass\vcat([-k*(u[2]-(200+273.15)); -fluid_props.rho(gasB)*fluid_props.enthalpy(gasB)*dV; -fluid_props.rho(gasB)*dV],
        (spec_MW[1:end-1]*1e-3).*[-r[1]+r[2];# production of CO in kg/s
            -r[2]-r[3];# production of CO2
            -2*r[1]-r[2]-3*r[3];# production of H2
            r[1]+r[3];# production of CH3OH
            r[2]+r[3]],# production of H2O,
        -fluid_props.Y(gasB)[spec_ind[1:end-1]].*fluid_props.rho(gasB)*dV)

    return nothing
end

function f2!(du, u, p, t)
    # f2 is the non-stiff portion
    Mass=setMassMat(u,p,t)
    gasA=p[1] # for passing to fluid props
    gasB=p[2] # for passing to fluid props
    V, dV = p[3](t) # function of t returning V and dV

    du .= Mass\vcat([-2.5e4*0*(u[2]-(25+273.15)); -fluid_props.rho(gasB)*fluid_props.enthalpy(gasB)*dV; -fluid_props.rho(gasB)*dV],
        zeros(5),
        -fluid_props.Y(gasB)[spec_ind[1:end-1]].*fluid_props.rho(gasB)*dV)

    return nothing
end


function f2_demo!(du, u, p, t)
    # f2 is the non-stiff portion
    Mass=setMassMat(u,p,t)
    gasA=p[1] # for passing to fluid props
    gasB=p[2] # for passing to fluid props
    V, dV = p[3](t) # function of t returning V and dV

    du .= Mass\vcat([-2.5e4*(u[2]-(25+273.15)); -fluid_props.rho(gasB)*fluid_props.enthalpy(gasB)*dV; -fluid_props.rho(gasB)*dV],
        zeros(5),
        -fluid_props.Y(gasB)[spec_ind[1:end-1]].*fluid_props.rho(gasB)*dV)

    return nothing
end


function f1!(du, u, p, t)
    # f1 is the stiff portions (ie reactions only)
    
    Mass=setMassMat(u,p,t)
    gasA=p[1] # for passing to fluid props
    gasB=p[2] # for passing to fluid props
    V, dV = p[3](t) # function of t returning V and dV
    rhocat=p[4]

    XA=fluid_props.X(gasA)[spec_ind]
    r= rhocat*[fluid_props.r1(u[2],u[1]/100e3,XA),
        fluid_props.r2(u[2],u[1]/100e3,XA), 
        fluid_props.r3(u[2],u[1]/100e3,XA)] # mol/s

    #tmp=Mass^-1
    #tmp[3,:] .*=-1

    du .= Mass\vcat([0;0;0],
        (spec_MW[1:end-1]*1e-3).*[-r[1]+r[2];# production of CO in kg/s
            -r[2]-r[3];# production of CO2
            -2*r[1]-r[2]-3*r[3];# production of H2
            r[1]+r[3];# production of CH3OH
            r[2]+r[3]],# production of H2O,
        zeros(Nspec-1))

    #println(tmp[2,4:8])
    #println(tmp[3,9:13])
    return nothing
end



function f1_test!(du, u, p, t)
    # f1 is the stiff portions (ie reactions only)
    
    Mass=setMassMat(u,p,t)
    gasA=p[1] # for passing to fluid props
    gasB=p[2] # for passing to fluid props
    V, dV = p[3](t) # function of t returning V and dV
    rhocat=p[4]

    XA=fluid_props.X(gasA)[spec_ind]
    r= rhocat*[fluid_props.r1(u[2],u[1]/100e3,XA),
        fluid_props.r2(u[2],u[1]/100e3,XA), 
        fluid_props.r3(u[2],u[1]/100e3,XA)] # mol/s


    dE=[intE0[4]*spec_MW[4]-intE0[1]*spec_MW[1]-2*intE0[3]*spec_MW[3];
        intE0[1]*spec_MW[1]+intE0[5]*spec_MW[5] - intE0[2]*spec_MW[2] - intE0[3]*spec_MW[3];
        intE0[4]*spec_MW[4] + intE0[5]*spec_MW[5] - intE0[2]*spec_MW[2] - 3* intE0[3]*spec_MW[3]]/1000; # all in units J/mol
    
    if t==0.0
        print(dE)
    end

    qrelease= sum(r .* dE) # in units J/s

    du .= Mass\vcat([-qrelease;0;0],
        (spec_MW[1:end-1]*1e-3).*[-r[1]+r[2];# production of CO in kg/s
            -r[2]-r[3];# production of CO2
            -2*r[1]-r[2]-3*r[3];# production of H2
            r[1]+r[3];# production of CH3OH
            r[2]+r[3]],# production of H2O,
        zeros(Nspec-1))

    #println(tmp[2,4:8])
    #println(tmp[3,9:13])
    return nothing
end

function τ_lookup(t)
    tfit_bnd=[5,20, 50,200,500,1000,2000]/1000 .+ 0.02
    τ=[  84.30031080810156
    148.23833599825355
    275.1187469982629
    337.2053964137515
    398.77163425418723
    592.4779359207819]/1000

    if t<tfit_bnd[1]
        return 0.1
    elseif t>tfit_bnd[end]
        return τ[end]
    else
        for i in eachindex(tfit_bnd[1:end-1])
            if t<tfit_bnd[i+1]
                return τ[i]
                break
            end
        end
    end

end

end