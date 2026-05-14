module thermochem_tools

    using PyCall
    using Printf
    using NLsolve
    const Ru=8.314 # J/mol K
    const meoh_prod=["CO","CO2","H2","H2O","N2","CH3OH"]

    TK=Vector{Float64}(undef, 0)
    KWGS=Vector{Float64}(undef, 0)
    KCOh=Vector{Float64}(undef, 0)

    ct=pyimport("cantera")

    include((@__DIR__)*"/pox_eqm.jl")
    include((@__DIR__)*"/meoh_eqm.jl")
    include((@__DIR__)*"/flammability.jl")
    include((@__DIR__)*"/Energy.jl")
    include((@__DIR__)*"/SRK_gas.jl")

    export ct

    function define_gas(mech::String)
        return ct.Solution(mech)
    end

    function comp_string(X::Vector{Float64},Specs::Vector{String})
        return join([@sprintf("%s:%.5f,", Specs[k], X[k]) for k in eachindex(Specs)])[1:end-1]
    end

    function set_state(gas, T,P,Comp::String)
        gas.TPX=(T,P,Comp)
        return nothing
    end





end