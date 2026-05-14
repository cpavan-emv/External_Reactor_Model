# D converts (rho, T, X) state -> (P, T, X) state
function Dcom!(D, gas)
    # IMPORTANT: assumes all inactive entries already zero
    D[1:2, 1:2] .= [gas.rho/gas.P  -gas.rho/gas.T;
                     0              1.0]
    D[3:end, 3:end] .= LinearAlgebra.diagm(ones(size(D, 1) - 2))
    for i in 1:size(D, 1)-2
        D[1, 2+i] = gas.rho / gas.MW_mix * (gas.MW_spec[i] - gas.MW_spec[end])
    end
    return nothing
end

function Dcom(gas::cantera_jl.gas_jl)
    N = gas.gas.Nspec + 1
    D = zeros(N, N)
    Dcom!(D, gas)
    return D
end

# B converts mole fractions to mass fractions: B * dX = dY
# Sum(Xi)=1 already imposed to eliminate Nth species
function Bcom!(B, gas)
    MW_mix  = gas.MW_mix
    spec_MW = gas.MW_spec
    X       = gas.X
    B .= 1/MW_mix * (LinearAlgebra.diagm(spec_MW[1:end-1]) -
         1/MW_mix * (transpose(spec_MW[1:end-1] .- spec_MW[end])) .* (spec_MW[1:end-1] .* X[1:end-1]))
    return nothing
end

function Bcom(gas)
    MW_mix  = gas.MW_mix
    spec_MW = gas.MW_spec
    X       = gas.X
    Nspec   = gas.gas.Nspec
    B = LinearAlgebra.diagm(ones(Nspec + 1))
    B[3:end, 3:end] .= 1/MW_mix * (LinearAlgebra.diagm(spec_MW[1:end-1]) -
        1/MW_mix * (transpose(spec_MW[1:end-1] .- spec_MW[end])) .* (spec_MW[1:end-1] .* X[1:end-1]))
    return B
end

# F converts d(rho T Y) to cv*dT + sum(yi * intE0_i)
function Fcom!(F, gas)
    h0 = gas.h0
    F .= LinearAlgebra.diagm(ones(size(F, 1)))
    F[2, 2] = gas.cv
    for i in 1:size(F, 1)-2
        F[2, 2+i] = h0[i] - h0[end]
    end
    return nothing
end

function Fcom(gas)
    Nspec = gas.gas.Nspec
    h0    = gas.h0
    F = LinearAlgebra.diagm(ones(Nspec + 1))
    F[2, 2] = gas.cv
    for i in 1:Nspec-1
        F[2, 2+i] = h0[i] - h0[end]
    end
    return F
end

function Fcom_vec(gas)
    Nspec = gas.gas.Nspec
    h0    = gas.h0
    F = zeros(Nspec + 1)
    F[2] = gas.cv
    for i in 1:Nspec-1
        F[2+i] = h0[i] - h0[end]
    end
    return F
end

# M is the mass matrix for a single region
function Mcom!(M, gas)
    # IMPORTANT: assumes all inactive entries already zero
    rho = gas.rho
    M[1:2, 1:2] .= [1.0 0.0; gas.int_nrg rho]
    M[3:end, 3:end] .= LinearAlgebra.diagm(ones(size(M, 1) - 2) * rho)
    M[3:end, 1] .= gas.Y[1:end-1]
    return nothing
end

function Mcom(gas)
    N = gas.gas.Nspec + 1
    M = zeros(N, N)
    Mcom!(M, gas)
    return M
end

# Couples clearance (B) and displaced (C) volumes under constant-pressure constraint
function Mass_BC(gasB, gasC, V_dV, scratch::ReactorScratch)
    Nr  = gasB.gas.Nspec + 1   # per-region state dimension [P, T, X1..X_{Nspec-1}]
    MCB = scratch.MCB
    FBD = scratch.FBD
    fill!(MCB, 0.0)
    fill!(FBD, 0.0)

    Mcom!((@view MCB[1:Nr, 1:Nr]), gasB)
    Mcom!((@view MCB[Nr+1:end, Nr+1:end]), gasC)

    # dV>0: flow B->C (B upstream); otherwise C->B
    gasup = V_dV[2] > 0 ? gasB : gasC
    Yup   = gasup.Y[1:end-1]

    # Modifications to region B equations
    MCB[1,    Nr+1]    += V_dV[1]
    MCB[2,    Nr+1]    += gasup.enthalpy * V_dV[1]
    MCB[3:Nr, Nr+1]   .+= Yup * V_dV[1]

    # Modifications to region C equations
    MCB[Nr+2,     Nr+1]    += -gasup.enthalpy
    MCB[Nr+3:end, Nr+1]   .+= -Yup

    FBD_B = Fcom(gasB) * Bcom(gasB) * Dcom(gasB)
    FBD_C = Fcom(gasC) * Bcom(gasC) * Dcom(gasC)
    FBD[1:Nr, 1:Nr]         .= FBD_B
    FBD[Nr+1:end, Nr+1:end] .= FBD_C

    mul!(scratch.BC_prod, MCB, FBD)

    # Replace region C continuity equation with constant-pressure constraint: P_B = P_C
    scratch.BC_prod[Nr+1, :]    .= 0.0
    scratch.BC_prod[Nr+1, 1]     = 1.0
    scratch.BC_prod[Nr+1, Nr+1] = -1.0

    return scratch.BC_prod
end

# Mass matrix for region A (always independent)
function Mass_A(gas)
    return Mcom(gas) * Fcom(gas) * Bcom(gas) * Dcom(gas)
end

# Full 3-region mass matrix: A independent, B and C coupled at constant pressure
function Mass_full(u, params::ReactorParams, t)
    scratch  = params.scratch
    gasses   = params.gasses
    V_dV     = params.Vfunc(t)
    N        = gasses[1].gas.Nspec + 1
    for i in 1:3
        offset   = N * (i - 1)
        X_partial = u[offset+3:offset+N]
        ct.setTPX(gasses[i], (u[offset+2], u[offset+1], [X_partial; 1.0 - sum(X_partial)]))
    end
    fill!(scratch.Mass, 0.0)
    scratch.Mass[1:N, 1:N]         .= Mass_A(gasses[1])
    scratch.Mass[N+1:end, N+1:end] .= Mass_BC(gasses[2], gasses[3], V_dV, scratch)
    return scratch.Mass
end
