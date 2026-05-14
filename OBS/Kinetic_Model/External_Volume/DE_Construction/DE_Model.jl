module DE_Model

# include the required dependencies
include((@__DIR__)*"/../../fluid_properties/fluid_props.jl")
include((@__DIR__)*"/Common_Matrices.jl")
include((@__DIR__)*"/Coupled_Chambers.jl")
include((@__DIR__)*"/CV_Reactor.jl")
include((@__DIR__)*"/Compression.jl")
include((@__DIR__)*"/Decoupled_Chambers.jl")
include((@__DIR__)*"/Intake_Exhaust.jl")

using LinearAlgebra, NLsolve
import OrdinaryDiffEq  as ODE

spec_ind=[0] # index of all species
spec_MW=[0.0] # molecular weight of all species
Nspec=0 # number of species
intE0=[0.0] # ref. state internal energy of each species (J/kg)
fluid_props=fluid_props
M_decoup2coup=Matrix{Float64}(undef,0,0)

function initialize_ideal_gas(ct_mech)
    gas=fluid_props.initialize_ideal_gas(ct_mech)
    if isnothing(fluid_props.kingas)
        fluid_props.kingas=gas
    end
    return gas
end

function set_gas_constants(gas, specs)
    # this is basically the initializer
    DE_Model.spec_ind=fluid_props.spec_inds(gas,specs) .+ 1
    DE_Model.spec_MW=fluid_props.MW_spec(gas)[spec_ind]
    DE_Model.Nspec=length(specs)
    DE_Model.M_decoup2coup=_decoupled2coupled_mat()
    # can use either h or u here since everything becomes deltas
    # stick with h because I like enthalpy of formation
    #DE_construct.intE0=fluid_props.u0(gas, specs) # ref. state internal energy of each species (J/kg)
    DE_Model.intE0=fluid_props.h0(gas, specs) # ref. state internal energy of each species (J/kg)
end



end
