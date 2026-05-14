# D converts (rho, T, X) state → (P, T, X) state
function Dcom!(D, gas, cfg::ModelConfig)
    # IMPORTANT: assumes all inactive entries already zero
    D[1:2, 1:2] .= [gas.rho/gas.P  -gas.rho/gas.T;
                     0              1.0]
    D[3:end, 3:end] .= LinearAlgebra.diagm(ones(size(D, 1) - 2))
    for i in 1:size(D, 1)-2
        D[1, 2+i] = gas.rho / gas.MW_mix * (cfg.spec_MW[i] - cfg.spec_MW[end])
    end
    return nothing
end

function Dcom(gas, cfg::ModelConfig)
    D = zeros(cfg.Nspec + 1, cfg.Nspec + 1)
    Dcom!(D, gas, cfg)
    return D
end

# B converts mole fractions to mass fractions: B * dX = dY
# Sum(Xi)=1 already imposed to eliminate Nth species
function Bcom!(B, gas, cfg::ModelConfig)
    MW = gas.MW_mix
    X  = gas.X[cfg.spec_ind]
    B .= 1/MW * (LinearAlgebra.diagm(cfg.spec_MW[1:end-1]) -
         1/MW * (transpose(cfg.spec_MW[1:end-1] .- cfg.spec_MW[end])) .* (cfg.spec_MW[1:end-1] .* X[1:end-1]))
    return nothing
end

function Bcom(gas, cfg::ModelConfig)
    MW = gas.MW_mix
    X  = gas.X[cfg.spec_ind]
    B  = LinearAlgebra.diagm(ones(cfg.Nspec + 1))
    B[3:end, 3:end] .= 1/MW * (LinearAlgebra.diagm(cfg.spec_MW[1:end-1]) -
        1/MW * (transpose(cfg.spec_MW[1:end-1] .- cfg.spec_MW[end])) .* (cfg.spec_MW[1:end-1] .* X[1:end-1]))
    return B
end

# F converts d(rho T Y) to cv*dT + sum(yi * intE0_i)
function Fcom!(F, gas, cfg::ModelConfig)
    F .= LinearAlgebra.diagm(ones(size(F, 1)))
    F[2, 2] = gas.cv
    for i in 1:size(F, 1)-2
        F[2, 2+i] = cfg.intE0[i] - cfg.intE0[end]
    end
    return nothing
end

function Fcom(gas, cfg::ModelConfig)
    F = LinearAlgebra.diagm(ones(cfg.Nspec + 1))
    F[2, 2] = gas.cv
    for i in 1:cfg.Nspec-1
        F[2, 2+i] = cfg.intE0[i] - cfg.intE0[end]
    end
    return F
end

function Fcom_vec(gas, cfg::ModelConfig)
    F = zeros(cfg.Nspec + 1)
    F[2] = gas.cv
    for i in 1:cfg.Nspec-1
        F[2+i] = cfg.intE0[i] - cfg.intE0[end]
    end
    return F
end

# M is the mass matrix for a single region
function Mcom!(M, gas, cfg::ModelConfig)
    # IMPORTANT: assumes all inactive entries already zero
    rho = gas.rho
    M[1:2, 1:2] .= [1.0 0.0; gas.int_nrg rho]
    M[3:end, 3:end] .= LinearAlgebra.diagm(ones(size(M, 1) - 2) * rho)
    M[3:end, 1] .= gas.Y[cfg.spec_ind[1:end-1]]
    return nothing
end

function Mcom(gas, cfg::ModelConfig)
    M = zeros(cfg.Nspec + 1, cfg.Nspec + 1)
    Mcom!(M, gas, cfg)
    return M
end

# Couples clearance (B) and displaced (C) volumes under constant-pressure constraint
function Mass_BC(gasB, gasC, V_dV, cfg::ModelConfig, scratch::ReactorScratch)
    Nr  = cfg.Nspec + 1   # per-region state dimension [P, T, X1..X_{Nspec-1}]
    MCB = scratch.MCB
    FBD = scratch.FBD
    fill!(MCB, 0.0)
    fill!(FBD, 0.0)

    Mcom!((@view MCB[1:Nr, 1:Nr]), gasB, cfg)
    Mcom!((@view MCB[Nr+1:end, Nr+1:end]), gasC, cfg)

    # dV>0: flow B→C (B upstream); otherwise C→B
    gasup = V_dV[2] > 0 ? gasB : gasC
    Yup   = gasup.Y[cfg.spec_ind[1:end-1]]

    # Modifications to region B equations
    MCB[1,    Nr+1]    += V_dV[1]
    MCB[2,    Nr+1]    += gasup.enthalpy * V_dV[1]
    MCB[3:Nr, Nr+1]   .+= Yup * V_dV[1]

    # Modifications to region C equations
    MCB[Nr+2,     Nr+1]    += -gasup.enthalpy
    MCB[Nr+3:end, Nr+1]   .+= -Yup

    FBD_B = Fcom(gasB, cfg) * Bcom(gasB, cfg) * Dcom(gasB, cfg)
    FBD_C = Fcom(gasC, cfg) * Bcom(gasC, cfg) * Dcom(gasC, cfg)
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
function Mass_A(gas, cfg::ModelConfig)
    return Mcom(gas, cfg) * Fcom(gas, cfg) * Bcom(gas, cfg) * Dcom(gas, cfg)
end

# Full 3-region mass matrix: A independent, B and C coupled at constant pressure
function Mass_full(u, params::ReactorParams, t)
    cfg     = params.config
    scratch = params.scratch
    gasses  = params.gasses
    V_dV    = params.Vfunc(t)
    N       = cfg.Nspec + 1
    for i in 1:3
        offset = N * (i - 1)
        ct.setTPX(gasses[i], (u[offset+2], u[offset+1], u[offset+3:offset+N]), cfg.spec_ind)
    end
    fill!(scratch.Mass, 0.0)
    scratch.Mass[1:N, 1:N]         .= Mass_A(gasses[1], cfg)
    scratch.Mass[N+1:end, N+1:end] .= Mass_BC(gasses[2], gasses[3], V_dV, cfg, scratch)
    return scratch.Mass
end
