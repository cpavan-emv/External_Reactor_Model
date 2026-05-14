module DE_Model

using LinearAlgebra, NLsolve
import OrdinaryDiffEq as ODE
import DifferentialEquations: terminate!

include((@__DIR__)*"/../../submodules/cantera_julia/src/cantera.jl")
const ct = cantera_jl

# Default ODE solver: Rodas5P with finite-difference Jacobians.
# AutoFiniteDiff is permanently required because the Cantera C-library interface
# is opaque to Julia's AD systems (ForwardDiff, Enzyme, Zygote) — same as PyCall was.
const _default_alg = ODE.Rodas5P(autodiff=ODE.AutoFiniteDiff())

# ── Mechanism constants (set once per mechanism, then passed through params) ──
struct ModelConfig
    spec_ind::Vector{Int}     # 1-based Cantera species indices for tracked species
    spec_MW::Vector{Float64}  # molecular weights (kg/kmol) for tracked species
    Nspec::Int                # number of tracked species (including bath gas)
    intE0::Vector{Float64}    # formation enthalpies h0 (J/kg) for tracked species
end

# ── Pre-allocated scratch buffers for the ODE hot path ───────────────────────
mutable struct ReactorScratch
    Mass::Matrix{Float64}    # 3N×3N  full mass matrix
    MCB::Matrix{Float64}     # 2N×2N  B+C block (before FBD transform)
    FBD::Matrix{Float64}     # 2N×2N  block-diagonal FBD transform for B+C
    BC_prod::Matrix{Float64} # 2N×2N  result of MCB*FBD (B+C mass matrix)
    RHS::Vector{Float64}     # 3N     right-hand side buffer
    r::Vector{Float64}       # 3      kinetic rates [r1, r2, r3]
end

function ReactorScratch(Nspec::Int)
    N = Nspec + 1
    ReactorScratch(zeros(3N, 3N), zeros(2N, 2N), zeros(2N, 2N), zeros(2N, 2N),
                   zeros(3N), zeros(3))
end

# ── Simulation parameters (one per simulation) ───────────────────────────────
struct ReactorParams{G, F<:Function}
    gasses::Vector{G}         # gas_jl objects, one per region (A, B, C)
    Vfunc::F                  # t → (V, dV)  displaced-volume kinematics
    rhocat::Float64           # catalyst effective density (kg/m³)
    T_walls::Vector{Float64}  # wall temperatures [T_A, T_B, T_C] (K)
    τ::Vector{Float64}        # heat-loss time constants [τ_A, τ_B, τ_C] (s)
    K_vals::Vector{Float64}   # pressure-drop coefficients:
                              #   [K_BA_rev, K_BA_fwd, K_IE_in, K_IE_out]
    config::ModelConfig
    scratch::ReactorScratch
end

function ReactorParams(gasses::Vector{G}, Vfunc::F, rhocat, T_walls, τ, K_vals,
                       config::ModelConfig) where {G, F<:Function}
    ReactorParams{G,F}(gasses, Vfunc, rhocat, T_walls, τ, K_vals, config,
                       ReactorScratch(config.Nspec))
end

# ── ODE mode hierarchy ────────────────────────────────────────────────────────
abstract type ReactorMode end

struct Decoupled <: ReactorMode end  # A and B+C volumes pressure-decoupled

struct Coupled <: ReactorMode        # A and B+C connected via finite pressure drop
    Vex::Float64                     # external reactor volume (same units as Vfunc output)
end

struct IntakeExhaust <: ReactorMode  # intake or exhaust valve open
    TPX_exterior::Tuple              # (T_K, P_Pa, X_vec_or_string) of the manifold
end

include((@__DIR__)*"/Common_Matrices.jl")
include((@__DIR__)*"/RHS_functions.jl")
include((@__DIR__)*"/Intake_Exhaust.jl")

# ── Module globals (set by set_gas_constants for legacy / interactive use) ────
spec_ind = Int[]
spec_MW  = Float64[]
Nspec    = 0
intE0    = Float64[]
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

# Returns a ModelConfig and also updates the legacy module globals.
function set_gas_constants(gas, specs)
    DE_Model.spec_ind = [ct.get_speciesIndex(gas.gas, n) + 1 for n in specs]
    DE_Model.spec_MW  = gas.MW_spec[DE_Model.spec_ind]
    DE_Model.Nspec    = length(specs)
    DE_Model.intE0    = gas.h0[DE_Model.spec_ind]
    return ModelConfig(DE_Model.spec_ind, DE_Model.spec_MW, DE_Model.Nspec, DE_Model.intE0)
end

function Vars_Eq_condition(u, _t, _integrator, inds)
    return u[inds[1]] - u[inds[2]]
end

end
