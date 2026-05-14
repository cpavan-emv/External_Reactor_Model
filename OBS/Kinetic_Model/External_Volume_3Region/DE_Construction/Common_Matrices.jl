# D is used to convert from (rho, T, X) to (P, T, X)
function Dcom!(D, gas)
    # IMPORTANT: This assumes all inactive entries are already initialized to zero
    D[1:2,1:2] .= [fluid_props.rho(gas)/fluid_props.P(gas) -fluid_props.rho(gas)/fluid_props.T(gas);
        0 1.0]
    D[3:end,3:end] .= LinearAlgebra.diagm(ones(size(D,1)-2))
    for i=1:size(D,1)-2
        D[1,2+i] = fluid_props.rho(gas)/fluid_props.MW_mix(gas)*(spec_MW[i]-spec_MW[end])
    end
    return nothing
end

function Dcom(gas)
    D=zeros(2+Nspec-1, 2+Nspec-1)
    Dcom!(D, gas)
    return D
end

# B is the coefficients to convert from mass fraction to mole fraction
# i.e. B * dX = dY
# Sum(Xi)=1 already imposed to eliminate Nth species
function Bcom!(B, gas)
    MW=fluid_props.MW_mix(gas)
    X=fluid_props.X(gas)[spec_ind]
    B .= 1/MW * (LinearAlgebra.diagm(spec_MW[1:end-1]) - 
        1/MW*(transpose(spec_MW[1:end-1] .- spec_MW[end])) .* (spec_MW[1:end-1].*X[1:end-1])) 
    return nothing
end

function Bcom(gas)
    MW=fluid_props.MW_mix(gas)
    X=fluid_props.X(gas)[spec_ind]
    B=LinearAlgebra.diagm(ones(length(X)+1))
    B[3:end,3:end] .=1/MW * (LinearAlgebra.diagm(spec_MW[1:end-1]) - 
        1/MW*(transpose(spec_MW[1:end-1] .- spec_MW[end])) .* (spec_MW[1:end-1].*X[1:end-1])) 
    return B
end

# F converts du to cv*dT + sum(yi * u0_i)
# use as F * d(rho T Y)  = d(rho e Y)  

function Fcom!(F, gas)
    F .= LinearAlgebra.diagm(ones(size(F,1)))
    F[2,2]=fluid_props.cv(gas)
    for i=1:size(F,1)-2
        F[2,2+i] = (intE0[i]-intE0[end])
    end
    return nothing
end

function Fcom(gas)
    F = LinearAlgebra.diagm(ones(2+(Nspec-1)))
    F[2,2]=fluid_props.cv(gas)
    for i=1:Nspec-1
        F[2,2+i] = (intE0[i]-intE0[end])
    end
    return F
end

function Fcom_vec(gas)
    # this version is the row that is non-identity to speed up computation
    F=zeros(2+(Nspec-1))
    F[2]=fluid_props.cv(gas)
    for i=1:Nspec-1
        F[2+i] = (intE0[i]-intE0[end])
    end
    return F
end
# M is the mass matrix of the original system for each region independently

function Mcom!(M, gas)
    # IMPORTANT: This assumes all inactive entries are already initialized to zero
    rho=fluid_props.rho(gas)
    M[1:2,1:2] .= [1.0 0.0; fluid_props.int_nrg(gas) rho]
    M[3:end,3:end] .= LinearAlgebra.diagm(ones(size(M,1)-2)*rho)
    M[3:end,1] .= fluid_props.Y(gas)[spec_ind[1:end-1]]
    return nothing
end

function Mcom(gas)
    M=zeros(2+Nspec-1, 2+Nspec-1)
    Mcom!(M, gas)
    return M
end

# this couples together the clearance and displaced volumes
# this is done by constant pressure assumption
function Mass_BC(gasB, gasC, V_dV)
    # initialize matrix
    N_2=2+Nspec-1
    N=N_2*2
    MCB=zeros(Float64,(N,N))
    # add default versions
    Mcom!((@view MCB[1:N_2,1:N_2]), gasB)
    Mcom!((@view MCB[N_2+1:end,N_2+1:end]), gasC)
    # for region B, sub in the flow equation
    # dV>0 means flow from B to C (B upstream), otherwise C to B (C upstream)
    V_dV[2]>0 ? gasup=gasB : gasup=gasC
    Yup=fluid_props.Y(gasup)[spec_ind[1:end-1]]
    
    # Modification to region B equations
    MCB[1,N_2+1] += V_dV[1] # modification to continuity
    MCB[2,N_2+1] += fluid_props.enthalpy(gasup)*V_dV[1] # modification to energy cons.
    MCB[3:N_2, N_2+1] .+= Yup*V_dV[1] # modification to energy cons.

    # Modification to region C equations
    MCB[N_2+2,N_2+1] += -fluid_props.enthalpy(gasup) # modification to energy cons.
    MCB[N_2+3:end, N_2+1] .+= -Yup # modification to energy cons.

    # add all the conversion matrices and multiply through
    FBD_B=Fcom(gasB)*Bcom(gasB)*Dcom(gasB)
    FBD_C=Fcom(gasC)*Bcom(gasC)*Dcom(gasC)
    FBD=zeros(Float64,(N,N))
    FBD[1:N_2,1:N_2] .= FBD_B
    FBD[N_2+1:end, N_2+1:end] .= FBD_C
    Mass=MCB*FBD

    # Now apply the fixed pressure requirement to replace region C continuity
    Mass[N_2+1,:] .= 0.0
    Mass[N_2+1,1] = 1.0
    Mass[N_2+1, N_2+1] = -1.0

    return Mass
end

# Mass matrix for region A is always the same
function Mass_A(gas)
    return Mcom(gas)*
        Fcom(gas)*
        Bcom(gas)*
        Dcom(gas)
end

# full mass matrix evolves A independently and couples B and C
# coupling of B and C is constant pressure
# Flow of mass/species is from C -> B during compression
# vice versa for expansion
function Mass_full(u,p,t)
    # p=([gasA, gasB, gasC], V_dVfunc)
    gasses=p[1]
    V_dV = p[2](t) # function of t returning V and dV
    N=1+Nspec
    # first set all the gas states
    for i=1:3
        offset=N*(i-1)
        fluid_props.setTPX(gasses[i], (u[offset+2], u[offset+1], u[offset+3:offset+N]), spec_ind)
    end
    Mass=zeros(Float64,(3*N, 3*N))
    Mass[1:N, 1:N] .= Mass_A(gasses[1]) 
    Mass[N+1:end, N+1:end] .= Mass_BC(gasses[2], gasses[3], V_dV)
    return Mass
end