# D is used to convert from (rho,u, T, X) to (P, u, T, X)
function Dcom!(D, gas)
    # IMPORTANT: This assumes all inactive entries are already initialized to zero
    D[1,1] = fluid_props.rho(gas)/fluid_props.P(gas)
    D[1,3] = -fluid_props.rho(gas)/fluid_props.T(gas);
    D[2:end,2:end] .= LinearAlgebra.diagm(ones(size(D,1)-1))
    for i=1:size(D,1)-3
        D[1,3+i] = fluid_props.rho(gas)/fluid_props.MW_mix(gas)*(spec_MW[i]-spec_MW[end])
    end
    return nothing
end

function Dcom(gas)
    D=zeros(2+Nspec, 2+Nspec)
    Dcom!(D, gas)
    return D
end

# B is the coefficients to convert from mass fraction to mole fraction
# i.e. B * dX = dY
# Sum(Xi)=1 already imposed to eliminate Nth species
function Bcom!(B, gas)
    MW=fluid_props.MW_mix(gas)
    X=fluid_props.X(gas)[spec_ind]
    B[4:end, 4:end] .= 1/MW * (LinearAlgebra.diagm(spec_MW[1:end-1]) - 
        1/MW*(transpose(spec_MW[1:end-1] .- spec_MW[end])) .* (spec_MW[1:end-1].*X[1:end-1])) 
    return nothing
end

function Bcom(gas)
    MW=fluid_props.MW_mix(gas)
    X=fluid_props.X(gas)[spec_ind]
    B=LinearAlgebra.diagm(ones(length(X)+2))
    B[4:end,4:end] .=1/MW * (LinearAlgebra.diagm(spec_MW[1:end-1]) - 
        1/MW*(transpose(spec_MW[1:end-1] .- spec_MW[end])) .* (spec_MW[1:end-1].*X[1:end-1])) 
    return B
end

# F converts du to cv*dT + sum(yi * u0_i)
# use as F * d(rho u T Y)  = d(rho u e Y)  
function Fcom!(F, gas)
    F .= LinearAlgebra.diagm(ones(size(F,1)))
    F[3,3]=fluid_props.cv(gas)
    for i=1:size(F,1)-3
        F[3,3+i] = (intE0[i]-intE0[end])
    end
    return nothing
end

function Fcom(gas)
    F = LinearAlgebra.diagm(ones(2+(Nspec)))
    F[3,3]=fluid_props.cv(gas)
    for i=1:size(F,1)-3
        F[3,3+i] = (intE0[i]-intE0[end])
    end
    return F
end

# M is the mass matrix of the original system for each region independently
function Mcom!(M, gas, u)
    # IMPORTANT: This assumes all inactive entries are already initialized to zero
    #rho=fluid_props.rho(gas)
    #M[:,1] .= [u; fluid_props.int_nrg(gas); fluid_props.Y(gas)[spec_ind[1:end-1]]]
    #M[:,2:end] .= LinearAlgebra.diagm(ones(size(M,1))*rho)
    rho=fluid_props.rho(gas)
    M[:,1] .= [1.0; fluid_props.int_nrg(gas); fluid_props.Y(gas)[spec_ind[1:end-1]]]
    M[:,2:end] .= LinearAlgebra.diagm(ones(size(M,1))*rho)
    M[1,2] =0.0
    return nothing
end

function Mcom(gas, u)
    M=zeros(2+Nspec-1, 3+Nspec-1)
    Mcom!(M, gas, u)
    return M
end

function Mass_com(gas, u)
    tmp=zeros(Float64,Nspec+2)
    tmp[1]=1.0
    return vcat(tmp', Mcom(gas, u)*
        Fcom(gas)*Bcom(gas)*Dcom(gas))
end

function D_hY(gas)
    D = LinearAlgebra.diagm(ones(Nspec))
    D[1,1]=-fluid_props.rho(gas)/fluid_props.T(gas);
    D[1,2:end] = fluid_props.rho(gas)/fluid_props.MW_mix(gas)*(spec_MW[1:end-1] .- spec_MW[end])
    return D
end

function D_hY!(D,gas)
    D[1,:] .= fluid_props.rho(gas) * [-1/fluid_props.T(gas);
        1 ./ fluid_props.MW_mix(gas)*(spec_MW[1:end-1] .- spec_MW[end])]
    return nothing
end

function D_hY_rowvec(gas)
    return fluid_props.rho(gas) * [-1/fluid_props.T(gas);
        1 ./ fluid_props.MW_mix(gas)*(spec_MW[1:end-1] .- spec_MW[end])]
end


function FB_hY(gas)
    FB = LinearAlgebra.diagm(ones(Nspec))
    FB[1,:] .= [fluid_props.cv(gas);(intE0[1:end-1] .- intE0[end])]
    MW=fluid_props.MW_mix(gas)
    X=fluid_props.X(gas)[spec_ind]
    FB[2:end, 2:end] .=1/MW * (LinearAlgebra.diagm(spec_MW[1:end-1]) - 
        1/MW*(transpose(spec_MW[1:end-1] .- spec_MW[end])) .* (spec_MW[1:end-1].*X[1:end-1])) 
    return FB
end
function FB_hY!(FB,gas)
    FB[1,:] .= [fluid_props.cv(gas);(intE0[1:end-1] .- intE0[end])]
    MW=fluid_props.MW_mix(gas)
    X=fluid_props.X(gas)[spec_ind]
    FB[2:end, 2:end] .=1/MW * (LinearAlgebra.diagm(spec_MW[1:end-1]) - 
        1/MW*(transpose(spec_MW[1:end-1] .- spec_MW[end])) .* (spec_MW[1:end-1].*X[1:end-1])) 
    return nothing
end