# ct must be defined in the calling script's scope before these functions are called
# (typically: const ct = DE_Model.ct)

function get_rhoU(gas, y, region)
    Nspec = gas.gas.Nspec
    N     = Nspec + 1
    rho_u = Matrix{Float64}(undef, size(y, 1), 2)
    for i in axes(y, 1)
        offset    = (region - 1) * N
        X_partial = y[i, offset+3:offset+N]
        ct.setTPX(gas, (y[i, offset+2], y[i, offset+1], [X_partial; 1.0 - sum(X_partial)]))
        rho_u[i, :] .= [gas.rho, gas.int_nrg]
    end
    return rho_u
end

function get_rhoU_single_vol(gas, y)
    Nspec = gas.gas.Nspec
    rho_u = Matrix{Float64}(undef, size(y, 1), 2)
    for i in axes(y, 1)
        X_partial = y[i, 3:Nspec+1]
        ct.setTPX(gas, (y[i, 2], y[i, 1], [X_partial; 1.0 - sum(X_partial)]))
        rho_u[i, :] .= [gas.rho, gas.int_nrg]
    end
    return rho_u
end

function Mtot(rhoA, rhoB, V)
    return V .* rhoB .+ rhoA
end

function Utot(rhoUA, rhoUB, V)
    return V .* prod(rhoUB, dims=2) .+ prod(rhoUA, dims=2)
end

function Ecomp(P, V)
    Ein = similar(P)
    Ein[1] = 0.0
    deltaE = 0.5 * (P[2:end] .+ P[1:end-1]) .* (V[2:end] .- V[1:end-1])
    for i in eachindex(Ein[1:end-1])
        Ein[i+1] = Ein[i] + deltaE[i]
    end
    Ein .*= -1
    return Ein
end
