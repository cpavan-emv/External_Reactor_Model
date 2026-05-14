module DE_Model

# include the required dependencies
include((@__DIR__)*"/../../fluid_properties/fluid_props.jl")
include((@__DIR__)*"/Conversion_Matrices.jl")
include((@__DIR__)*"/RHS_functions.jl")
include((@__DIR__)*"/Spatial_Discretization.jl")
include((@__DIR__)*"/SS_Model.jl")
include((@__DIR__)*"/Trans_Model.jl")


using LinearAlgebra, NLsolve
import OrdinaryDiffEq  as ODE
import SparseArrays

spec_ind=[0] # index of all species
spec_MW=[0.0] # molecular weight of all species
Nspec=0 # number of species
intE0=[0.0] # ref. state internal energy of each species (J/kg)
fluid_props=fluid_props
xb=[0.0]
MapMat=[]
NvNc=[]


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
    #DE_Model.M_decoup2coup=_decoupled2coupled_mat()
    # can use either h or u here since everything becomes deltas
    # stick with h because I like enthalpy of formation
    #DE_construct.intE0=fluid_props.u0(gas, specs) # ref. state internal energy of each species (J/kg)
    DE_Model.intE0=fluid_props.h0(gas, specs) # ref. state internal energy of each species (J/kg)
end

function Vars_Eq_condition(u,t,integrator, inds)
    return u[inds[1]] - u[inds[2]]
end

function setup_grid(lims, Ncell)
    # this sets up linear spacing
    DE_Model.xb=collect(range(lims[1],lims[2], Ncell+1))
    return nothing
end

function setup_grid(xbou)
    # this allows custom spacing
    DE_Model.xb=xbou
    return nothing
end

function x_dx()
    dx=xb[2:end]-xb[1:end-1]
    x= 0.5*(xb[2:end]+xb[1:end-1])
    return x, dx
end

# define 2 orderings of parameters:
# first is "by Cell" where we put all variables for cell 1, then all for 2, etc.
# second is "by Var" where we put all values for variable 1, then all for 2, etc.
function MapMat_Cell2Var(Nvar, Ncell)
    # do this using a mapping matrix, which will be sparse
    # construct by R,C, V format
    N=Nvar*Ncell
    R=collect(1:N)
    C=similar(R)
    V=ones(Float64, N)
    for i=1:N
        # should re-order as 1,Nvar+1, 2*Nvar+1, etc.
        # then 2, Nvar+2, 2*Nvar+2, etc.
        C[i]=(floor(Int, (i-1)/Ncell)+1)+ # this is the index of the variable (1,2,3,4,...)+
            (mod(i-1,Ncell)*Nvar) # and this is the cell number to pull it from (multiplied by the length of each cell)
    end
    DE_Model.MapMat=SparseArrays.sparse(R, C, V)
    DE_Model.NvNc=[Nvar,Ncell]
end

function Cell2Var(state)
    # N is tupe of Nvar, Ncell
    return DE_Model.MapMat*state
end

function Var2Cell(state)
    # N is tupe of Nvar, Ncell
    return DE_Model.MapMat'*state
end




end