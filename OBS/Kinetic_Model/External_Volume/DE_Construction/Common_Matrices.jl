function Dcom!(D, gas)
    # IMPORTANT: This assumes all inactieve entries are already initialized to zero
    D[1:2,1:2] .= [fluid_props.rho(gas)/fluid_props.P(gas) -fluid_props.rho(gas)/fluid_props.T(gas);
        0 1.0]
    D[3:end,3:end] .= LinearAlgebra.diagm(ones(size(D,1)-2))
    for i=1:size(D,1)-2
        D[1,2+i] = fluid_props.rho(gas)/fluid_props.MW_mix(gas)*(spec_MW[i]-spec_MW[end])
    end
    return nothing
end
function Dcom(gas)
    # D is used to convert from (rho, T, X) to (P, T, X)
    D=zeros(2+Nspec-1, 2+Nspec-1)
    Dcom!(D, gas)
    return D
end

function Bcom!(B, gas)
    # B is the coefficients to convert from mass fraction to mole fraction
    MW=fluid_props.MW_mix(gas)
    X=fluid_props.X(gas)[spec_ind]
    B .= 1/MW * (LinearAlgebra.diagm(spec_MW[1:end-1]) - 
        1/MW*(transpose(spec_MW[1:end-1] .- spec_MW[end])) .* (spec_MW[1:end-1].*X[1:end-1])) 
    return nothing
end

function Bcom(gas)
    # B is the coefficients to convert from mass fraction to mole fraction
    # i.e. B * dX = dY
    MW=fluid_props.MW_mix(gas)
    X=fluid_props.X(gas)[spec_ind]
    return 1/MW * (LinearAlgebra.diagm(spec_MW[1:end-1]) - 
        1/MW*(transpose(spec_MW[1:end-1] .- spec_MW[end])) .* (spec_MW[1:end-1].*X[1:end-1])) 
end

function Fcom!(F, gas)
    # F converts du to cv*dT + sum(yi * u0_i)
    # use as F * d(rho T Y)  = d(rho e Y)  
    F .= LinearAlgebra.diagm(ones(size(F,1)))
    F[2,2]=fluid_props.cv(gas)
    for i=1:size(F,1)-2
        F[2,2+i] = (intE0[i]-intE0[end])
    end
    return nothing
end

function Fcom(gas)
    # F converts du to cv*dt + sum(yi * u0_i)
    F = LinearAlgebra.diagm(ones(2+(Nspec-1)))
    F[2,2]=fluid_props.cv(gas)
    for i=1:Nspec-1
        F[2,2+i] = (intE0[i]-intE0[end])
    end
    return F
end

function Fcom_vec(gas)
    # F converts du to cv*dt + sum(yi * u0_i)
    # this version is the row that is non-identity to speed up computation
    F=zeros(2+(Nspec-1))
    F[2]=fluid_props.cv(gas)
    for i=1:Nspec-1
        F[2+i] = (intE0[i]-intE0[end])
    end
    return F
end

function _decoupled2coupled_mat()
    M=zeros(Nspec*2+1,Nspec*2+2)
    M[1:2,1:2] .= [1 0; 0 1]
    M[3,Nspec+3] = 1
    M[4:4+Nspec-2,3:3+Nspec-2] .=LinearAlgebra.diagm(ones(Nspec-1))
    M[4+Nspec-1:end, 3+Nspec+1:end] .=LinearAlgebra.diagm(ones(Nspec-1))
    return M
end

function DC2CC(u)
    return M_decoup2coup*u
end

function CC2DC(u)
    uout=M_decoup2coup'*u
    uout[Int(length(uout)/2)+1]=uout[1]
    return uout
end