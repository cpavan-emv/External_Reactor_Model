function get_rhoU(gas,y, region)
    rho_u=Matrix{Float64}(undef, size(y,1),2)
    for i in axes(y,1)
        fluid_props.setTPX(gas, 
            (y[i,1+region], y[i,1], y[i,4+(region-1)*(DE_construct.Nspec-1):4+region*(DE_construct.Nspec-1)-1]), DE_construct.spec_ind)

        rho_u[i,:] .= [fluid_props.rho(gas), fluid_props.int_nrg(gas)]
    end
    return rho_u
end

function Mtot(rhoA, rhoB, V)
    return V.*rhoB .+ rhoA
end

function Utot(rhoUA, rhoUB, V)
    return V.*prod(rhoUB, dims=2) + prod(rhoUA, dims=2)
end

function Ecomp(P, V)
    Ein=similar(P)
    Ein[1]=0.0
    deltaE = 0.5*(P[2:end]+P[1:end-1]) .* (V[2:end]-V[1:end-1])
    for i in eachindex(Ein[1:end-1])
        Ein[i+1] =Ein[i]+deltaE[i]
    end
    Ein*=-1
    return Ein
end