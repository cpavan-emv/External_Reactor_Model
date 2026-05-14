function G(gas, T,P,spec)
    gas.TPX=T,P, @sprintf("%s:1", spec)
    return gas.enthalpy_mole - gas.entropy_mole * T
end

# equilibrium constants
function K_CO2hyd(gas, T,P)
    dG=G(gas, T,P,"CH3OH")+G(gas, T,P,"H2O")-
        G(gas, T,P,"CO2")-3*G(gas, T,P,"H2")
    return exp(-dG/(Ru*1e3*T))*(100e3/P)^2 # bar^-2
end
function K_COhyd(gas, T,P)
    dG=G(gas, T,P,"CH3OH")-
        G(gas, T,P,"CO")-2*G(gas, T,P,"H2")
    return exp(-dG/(Ru*1e3*T))*(100e3/P)^2 # bar^-2
end
function K_WGS(gas, T,P)
    dG=G(gas, T,P,"CO2")+G(gas, T,P,"H2") - 
    G(gas, T,P,"CO")-G(gas, T,P,"H2O")
    return exp(-dG/(Ru*1e3*T))
end
# fits from Bennekom paper

# function K_COhyd(gas, T,P)
#     print("Here Once Again\n")
#     return 10^((5139 ./ T .- 12.621))#*(100e3/P)^2 # bar^-2
# end
# function K_WGS(gas, T,P)
#     return 10^(-(-2073 ./ T .+ 2.029))
# end

function quad_interp(x,y, xq)
    if any(xq .== x)
        return last(y[xq.==x])
    end
    Nx=length(x)
    inds=1:Nx
    indclose=inds[argmin(abs.(xq .- x))]
    rng=Vector{Int64}(undef,3)
    if indclose==1
        rng .= [1,2,3]
    elseif indclose==Nx
        rng .= [Nx-2, Nx-1, Nx]
    else
        rng .= [indclose-1,indclose, indclose+1]
    end

    C=[x[rng].^2 x[rng] ones(3)] \ y[rng]
    return C[1]*xq^2+C[2]*xq+C[3]
end

function count_atoms!(nA,x)
    #=
    assumes order is [CO, CO2, H2, H2O, N2, MeOH]
    same as the constant "meoh_prod"
    returns [C,H,O,N] (in-place)
    =#
    nA[1:4]=[x[1]+x[2]+x[6]; # Carbon
        2*(x[3]+x[4])+4*x[6]; # hydrogen
        x[1]+2*x[2]+x[4]+x[6]; # oxygen
        2*x[5]] # nitrogen
    return nothing
end

function count_atoms(x)
    nA=Vector{Float64}(undef,4)
    count_atoms!(nA,x)
    return nA
end

function MeOH_eqm_SP_res!(F,u,p, real_gas::Bool)
    # u is the vector of molar densities
    # p (parameters) is [T, P, xfresh, gas_object]
    T=p[1]; P=p[2]
    nA0=Vector{Float64}(undef,4)
    count_atoms!(nA0,p[3])
    x=u/sum(u)
    if real_gas
        tmp=fugacity_coeff(P,T,construct_xvec(x,meoh_prod))
        fug_coeff=[tmp[SRK_spec[s]] for s in meoh_prod]
        x .*= fug_coeff
        Kw=quad_interp(thermochem_tools.TK, thermochem_tools.KWGS, T)
        Kco=quad_interp(thermochem_tools.TK, thermochem_tools.KCOh, T)
    else
        Kw=(K_WGS(p[4],T,100e3))
        Kco=K_COhyd(p[4],T,100e3)
    end

    # 3 equilibrium reactions - only need 2 of them
    #F[4]= x[4]*x[6] - (x[2]*x[3]^3)*K_COhyd(T,100e3)*(P)^2
    F[5]= x[6] - (x[1]*x[3]^2)*Kco*(P/1e5)^2
    F[6]= x[2]*x[3] - (x[1]*x[4])*Kw
    # 4 atomic balance
    F[1:4] .= [u[1] + u[2] + u[6] - nA0[1]; # C
        2*u[3] + 2*u[4] + 4*u[6] - nA0[2]; # H
        u[1] + 2*u[2] + u[4] + u[6] - nA0[3]; # O
        2*u[5] - nA0[4]] # N2
end

function meoh_eq_SP(gas, Xin, Treact, Preact, real_gas=false)
    uin=copy(Xin)
    uin[4]=0.0 # remove liquid water
    uin[6]=0.0 # remove any methanol
    uin/=sum(uin) # normalize to 1 mole inlet
    comp=comp_string(uin, meoh_prod)
    gas.TPX=(Treact,Preact,comp)
    MW0=gas.mean_molecular_weight
    param=(Treact,Preact, uin, gas)
    fun(F,u)=MeOH_eqm_SP_res!(F, u, param, real_gas)
    sol=nlsolve(fun, uin)
    if any(sol.zero .< 0.0)
        # try again with different IC
        uin=2.0 / 6 *ones(6)
        sol=nlsolve(fun, uin)
    end
    comp=comp_string(sol.zero, meoh_prod)
    gas.TPX=(Treact,Preact,comp)
    MW1=gas.mean_molecular_weight

    return sol.zero/sum(sol.zero), MW0/MW1
end

function MeOH_eqm_recirc_res!(F,u,p, real_gas)
    # u is the vector of molar densities plus mfr 
    # p (parameters) is [T, P, xfresh, frac_purge, gas_object]
    T=p[1]; P=p[2]

    uliq=u .* [0,0,0,1.0,0,1.0] # liquid produced
    ugas=u-uliq # gas produced
    upurge=ugas * p[4] # purged gas
    urecirc = ugas-upurge # reciculated gas
    uin=urecirc + p[3] # gas into reactor

    nA0=Vector{Float64}(undef,4)
    count_atoms!(nA0,uin) # atom count into reactor
    nA=similar(nA0)
    count_atoms!(nA,u) # atom count out of reactor
    x=u/sum(u) # mole fractions in outlet
    # 4 atomic balance

    if real_gas
        tmp=fugacity_coeff(P,T,construct_xvec(x,meoh_prod))
        fug_coeff=[tmp[SRK_spec[s]] for s in meoh_prod]
        x .*= fug_coeff
        Kw=quad_interp(thermochem_tools.TK, thermochem_tools.KWGS, T)
        Kco=quad_interp(thermochem_tools.TK, thermochem_tools.KCOh, T)
    else
        Kw=(K_WGS(p[4],T,100e3))
        Kco=K_COhyd(p[4],T,100e3)
    end

    F[1:4] .= nA0-nA # N2
    # 2 equilibrium reactions
    F[5]= x[6] - (x[1]*x[3]^2)*Kco*(P/1e5)^2
    F[6]= x[2]*x[3] - (x[1]*x[4])*Kw
end

function meoh_eq_recirc(gas, Xin, Treact, Preact, fpurge, u0=[0.0], real_gas=false)
    uin=copy(Xin)
    uin[4]=0.0 # remove liquid water
    uin[6]=0.0 # remove any methanol
    uin/=sum(uin) # normalize to 1 mole inlet
    comp=comp_string(uin, meoh_prod)
    gas.TPX=(Treact,Preact,comp)
    MW0=gas.mean_molecular_weight
    param=(Treact,Preact, uin, fpurge, gas)
    fun(F,u)=MeOH_eqm_recirc_res!(F, u, param, real_gas)
    usol=similar(uin)
    if all(u0 .== 0.0)
        u0 .= uin
    end
    sol=nlsolve(fun, u0)
    usol .= sol.zero
    uout=([0;0;0;1;0;1]+fpurge*[1;1;1;0;1;0]) .* usol
    comp=comp_string(uout, meoh_prod)
    gas.TPX=(Treact,Preact,comp)
    MW1=gas.mean_molecular_weight

    return usol/sum(usol), MW0/MW1
end