##############################
#=
The functions in this file are for describing the case of compression with 2 volumes
Volume A is the volume of fixed size and contains the catalyst (reactions allowed)
Volume B is the volume of variable size and contains no catalyst (no reactions)
The two volumes are constrained to have the same pressure with mass/energy flow between them
=#
##############################


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

# function Dmat(gasA, gasB)
#     # D is the coefficients to convert rho to a function of P and T
#     return [fluid_props.rho(gasA)/fluid_props.P(gasA) -fluid_props.rho(gasA)/fluid_props.T(gasA)    0;
#         fluid_props.rho(gasB)/fluid_props.P(gasB)   0   -fluid_props.rho(gasB)/fluid_props.T(gasB);
#         0   1.0    0;
#         0   0   1.0]
# end

# function Bmat(gas)
#     # B is the coefficients to convert from mass fraction to mole fraction
#     MW=fluid_props.MW_mix(gas)
#     X=fluid_props.X(gas)[spec_ind]
#     return 1/MW * (LinearAlgebra.diagm(spec_MW[1:end-1]) - 
#         1/MW*(transpose(spec_MW[1:end-1] .- spec_MW[end])) .* (spec_MW[1:end-1].*X[1:end-1])) 
# end

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

# function Dfull(gasA, gasB)
#     # D converts rho to P and T
#     D=zeros(4+2*(Nspec-1), 3+2*(Nspec-1))
#     D[1:4, 1:3] .= Dmat(gasA, gasB)
#     D[5:end,4:end] .= LinearAlgebra.diagm(ones(2*(Nspec-1)))
#     for i=1:Nspec-1
#         D[1,3+i] = fluid_props.rho(gasA)/fluid_props.MW_mix(gasA)*(spec_MW[i]-spec_MW[end])
#         D[2,3+Nspec-1+i] = fluid_props.rho(gasB)/fluid_props.MW_mix(gasB)*(spec_MW[i]-spec_MW[end])
#     end
#     return D
# end

function D_CC(gasA, gasB)
    DA=Dcom(gasA)
    DB=Dcom(gasB)
    D=zeros(4+2*(Nspec-1), 3+2*(Nspec-1))
    # for region A
    D[1,1]=DA[1,1] # pressure coefficient
    D[1,2]=DA[1,2] # region A temperature coefficient
    D[1,4:4+Nspec-2]=DA[1,3:end] # region A species coefficients
    D[3,2]=DA[2,2] # region A temperature coefficient
    D[5:5+Nspec-2,4:4+Nspec-2]=DB[3:end,3:end] # should be a diagonal identity matrix
     # for region B
    D[2,1]=DB[1,1] # pressure coefficient
    D[2,3]=DB[1,2] # region B temperature coefficient
    D[2,:4+Nspec-1:end]=DB[1,3:end]
    D[4,3]=DA[2,2] # region A temperature coefficient
    D[5+Nspec-1:end,4+Nspec-1:end]=DB[3:end,3:end] # should be a diagonal identity matrix
    return D
end


# function Bfull(gasA, gasB)
#     # B converts mass fraction to mole fraction
#     B=LinearAlgebra.diagm(ones(4+2*(Nspec-1)))
#     B[5:5+(Nspec-2),5:5+(Nspec-2)] .= Bmat(gasA)
#     B[5+Nspec-1:end,5+Nspec-1:end] .= Bmat(gasB)
#     return B
# end

function B_CC(gasA, gasB)
    # B converts mass fraction to mole fraction
    B=LinearAlgebra.diagm(ones(4+2*(Nspec-1)))
    B[5:5+(Nspec-2),5:5+(Nspec-2)] .= Bcom(gasA)
    B[5+Nspec-1:end,5+Nspec-1:end] .= Bcom(gasB)
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

# function F_CC(gasA, gasB)
#     F=zeros(4+2*(Nspec-1), 4+2*(Nspec-1))
#     F[[1; 3; collect(5:5+Nspec-2)], [1;3;collect(5:5+Nspec-2)]] .= Fcom(gasA)
#     F[[2; 4; collect(5+Nspec-1:size(F,1))], [2; 4; collect(5+Nspec-1:size(F,1))]] .= Fcom(gasB)
#     return F
# end

function F_CC(gasA, gasB)
    F=LinearAlgebra.diagm(ones(4+2*(Nspec-1)))
    F[3,[1;3;collect(5:5+Nspec-2)]] .= Fcom_vec(gasA)
    F[4,[2; 4; collect(5+Nspec-1:size(F,1))]] .= Fcom_vec(gasB)
    return F
end


function setMassMat_CC(u,p,t)
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
         F_CC(gasA,gasB)*
         B_CC(gasA,gasB)*
         D_CC(gasA,gasB)
end

function f_CC!(du, u, p, t)
    # this is the entire RHS of the function
    # should only be used with stiff equation solvers
    Mass=setMassMat_CC(u,p,t)
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
    r[isnan.(r)] .= 0.0 # over-ride for case of pure N2 causing NANs

    τ=50/1000#τ_lookup(t)
    k=fluid_props.rho(gasA)*fluid_props.cp(gasA)/τ
    #println(k)


    du .= Mass\vcat([-k*(u[2]-(Tcat)); -fluid_props.rho(gasB)*fluid_props.enthalpy(gasB)*dV; -fluid_props.rho(gasB)*dV],
        (spec_MW[1:end-1]*1e-3).*[-r[1]+r[2];# production of CO in kg/s
            -r[2]-r[3];# production of CO2
            -2*r[1]-r[2]-3*r[3];# production of H2
            r[1]+r[3];# production of CH3OH
            r[2]+r[3]],# production of H2O,
        -fluid_props.Y(gasB)[spec_ind[1:end-1]].*fluid_props.rho(gasB)*dV)

    return nothing
end

function evolve_CC(u0,param,tsave, condition=nothing)
    prob=ODE.ODEProblem(DE_Model.f_CC!,u0,(tsave[1],tsave[end]))
    if !isnothing(condition)
        cb=ODE.ContinuousCallback(condition, (integrator)->ODE.terminate!(integrator))
        sol = ODE.solve(prob,ODE.RadauIIA5(autodiff=false),
            p=param, reltol=1e-8, abstol=1e-10, saveat=tsave, callback=cb)
    else
        sol = ODE.solve(prob,ODE.RadauIIA5(autodiff=false),
            p=param, reltol=1e-8, abstol=1e-10, saveat=tsave)
    end
    y=permutedims(stack(sol.u))
    return sol.t, y
end