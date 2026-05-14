##############################
#=
The functions in this file are for describing the case of a constant volume reactor
The system has no inlets or outlets and reactions are allowed
System is of form:
M ∂_t [rho e Y] = RHS
=#
##############################

function M_CV(gasA)
    # M is the coefficient matrix of the base equation
    M=zeros(2+(Nspec-1),2+(Nspec-1))
    Y=fluid_props.Y(gasA)[spec_ind]
    rho=fluid_props.rho(gasA)
    M[1:2,1:2] .= [1.0 0.0; fluid_props.int_nrg(gasA) rho]
    for i=1:Nspec-1
        M[2+i, 2+i] = rho
        M[2+i, 1] = Y[i]
    end
    return M
end

# Matrix D-F are all the same as the "common" version
function D_CV(gasA)
    return Dcom(gasA)
end

function B_CV(gasA)
    B=LinearAlgebra.diagm(ones(2+Nspec-1))
    B[3:end, 3:end].=Bcom(gasA)
    return B
end

function F_CV(gasA)
    return Fcom(gasA)
end


function setMassMat_CV(u,p,t)
    gasA=p[1] # for passing to fluid props
    fluid_props.setTPX(gasA, (u[2], u[1], u[3:end]), spec_ind)
    return M_CV(gasA)*
         F_CV(gasA)*
         B_CV(gasA)*
         D_CV(gasA)
end

function f_CV!(du, u, p, t)
    # this is the entire RHS of the function
    # should only be used with stiff equation solvers
    Mass=setMassMat_CV(u,p,t)
    gasA=p[1] # for passing to fluid props
    rhocat=p[4]
    Tcat=p[5]

    XA=fluid_props.X(p[1])[spec_ind]
    # r= rhocat*[fluid_props.r1(u[2],u[1]/100e3,XA),
    #     fluid_props.r2(u[2],u[1]/100e3,XA), 
    #     fluid_props.r3(u[2],u[1]/100e3,XA)] # mol/s
    # # fix catalyst temp
    r= rhocat*[fluid_props.r1(Tcat,u[1]/100e3,XA),
        fluid_props.r2(Tcat,u[1]/100e3,XA), 
        fluid_props.r3(Tcat,u[1]/100e3,XA)] # mol/s
    r[isnan.(r)] .= 0.0 # over-ride for case of pure N2 causing NANs

    τ=50/1000
    k=fluid_props.rho(gasA)*fluid_props.cp(gasA)/τ

    du .= Mass\vcat([0;-k*(u[2]-(Tcat))],
        (spec_MW[1:end-1]*1e-3).*[-r[1]+r[2];# production of CO in kg/s
            -r[2]-r[3];# production of CO2
            -2*r[1]-r[2]-3*r[3];# production of H2
            r[1]+r[3];# production of CH3OH
            r[2]+r[3]],# production of H2O,
            )

    return nothing
end







# function τ_lookup(t)
#     tfit_bnd=[5,20, 50,200,500,1000,2000]/1000 .+ 0.02
#     τ=[  84.30031080810156
#     148.23833599825355
#     275.1187469982629
#     337.2053964137515
#     398.77163425418723
#     592.4779359207819]/1000

#     if t<tfit_bnd[1]
#         return 0.1
#     elseif t>tfit_bnd[end]
#         return τ[end]
#     else
#         for i in eachindex(tfit_bnd[1:end-1])
#             if t<tfit_bnd[i+1]
#                 return τ[i]
#                 break
#             end
#         end
#     end

# end
