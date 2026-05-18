module DE_Model

using LinearAlgebra, NLsolve
import OrdinaryDiffEq as ODE
import OrdinaryDiffEqFIRK: RadauIIA5
import DifferentialEquations: terminate!

include((@__DIR__)*"/../../submodules/cantera_julia/src/cantera.jl")
const ct = cantera_jl

# Default ODE solver: Rodas5P with finite-difference Jacobians.
# AutoFiniteDiff is permanently required because the Cantera C-library interface
# is opaque to Julia's AD systems (ForwardDiff, Enzyme, Zygote) -- same as PyCall was.
#const _default_alg = ODE.Rodas5P(autodiff=ODE.AutoFiniteDiff())
const _default_alg = RadauIIA5(autodiff=ODE.AutoFiniteDiff())

# Pre-allocated scratch buffers for the ODE hot path
mutable struct ReactorScratch
    Mass::Matrix{Float64}    # 2N x 2N  full mass matrix
    RHS::Vector{Float64}     # 2N     right-hand side buffer
    r::Vector{Float64}       # 3      kinetic rates [r1, r2, r3]
end

function ReactorScratch(Nspec::Integer)
    N = Nspec + 1
    ReactorScratch(zeros(2N, 2N),zeros(2N), zeros(3))
end

# Simulation parameters (one per simulation).
struct ReactorParams{G, F<:Function}
    gases::Vector{G}         # gas_jl objects, one per region (A, B, C)
    Vfunc::F                  # t -> (V, dV)  displaced-volume kinematics (non-dim by displacement volumee)
    beta::Float64             # ratio of exgternal volume to displacement volume
    Vd::Float64               # displacement volume
    rhocat::Float64           # catalyst effective density (kg/m3)
    T_walls::Vector{Float64}  # wall temperatures [T_A, T_B, T_C] (K)
    tau::Vector{Float64}      # heat-loss time constants [tau_A, tau_B, tau_C] (s)
    Cv_vals::Vector{Float64}  # valve Cv [gal/min/√psi]: [Cv_BA_comp, Cv_BA_exp, Cv_IE_exh, Cv_IE_int]
    scratch::ReactorScratch
end

function ReactorParams(gases::Vector{G}, Vfunc::F, beta, Vd, rhocat, T_walls, tau, Cv_vals
                       ) where {G, F<:Function}
    ReactorParams{G,F}(gases, Vfunc, beta, Vd, rhocat, T_walls, tau, Cv_vals,
                       ReactorScratch(gases[1].gas.Nspec))
end

# ODE mode hierarchy
abstract type ReactorMode end

struct Decoupled <: ReactorMode end  # Separate external volume and cylinder

struct Coupled <: ReactorMode end      # External volume and cylinder connected                  

struct IntakeExhaust <: ReactorMode  # intake or exhaust valve open
    TPX_exterior::Tuple              # (T_K, P_Pa, X_vec_or_string) of the manifold
end

include((@__DIR__)*"/Common_Matrices.jl")
include((@__DIR__)*"/RHS_functions.jl")
include((@__DIR__)*"/Intake_Exhaust.jl")

# Module globals
_kinetics_initialized = false

function initialize_ideal_gas(ct_mech)
    ct.load_MeOH_kinetics()
    if !DE_Model._kinetics_initialized
        # load_MeOH_kinetics uses Base.include, creating bindings in a newer world.
        # invokelatest lets us call those bindings from here.
        Base.invokelatest() do
            ct.MeOH_kinetics.initialize_MeOH_kinetics(ct_mech)
        end
        DE_Model._kinetics_initialized = true
    end
    return ct.initialize_gas_jl(ct_mech)
end

# Used as a callback function
# triggers when the two variables are equal
# ex. when pressure in cylinder and external volume match
function Vars_Eq_condition(u, _t, _integrator, inds)
    return u[inds[1]] - u[inds[2]]
end

end
