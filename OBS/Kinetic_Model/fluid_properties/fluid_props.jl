module fluid_props

    using PyCall
    using Printf
    include((@__DIR__)*"/kinetics.jl")

    ct=nothing
    kingas=nothing

    function initialize_ideal_gas(ct_mech)
        fluid_props.ct=pyimport("cantera")
        return fluid_props.ct.Solution(ct_mech)
    end

    function initialize_kinetcs(ct_mech)
        fluid_props.ct=pyimport("cantera")
        return fluid_props.ct.Solution(ct_mech)
    end

    # ideal gas setter
    function setTPX(gas::PyObject, TPX::Tuple)
        gas.TPX=TPX
    end
    function setHPY(gas::PyObject, HPY::Tuple)
        gas.HPY=HPY
    end
    function setHPX(gas::PyObject, HPX::Tuple)
        gas.HPX=HPX
    end

    function setTPX(gas::PyObject,TPX::Tuple, spec_ind::Vector{Int})
        X=zeros(length(gas.X))
        X[spec_ind[1:end-1]].=TPX[end]
        X[spec_ind[end]] = 1-sum(X)
        setTPX(gas,(TPX[1],TPX[2],X))
        return nothing
    end

    function setHPY(gas::PyObject,HPY::Tuple, spec_ind::Vector{Int})
        Y=zeros(length(gas.Y))
        Y[spec_ind[1:end-1]].=HPY[end]
        Y[spec_ind[end]] = 1-sum(Y)
        setHPY(gas,(HPY[1],HPY[2],Y))
        return nothing
    end

    function setHPX(gas::PyObject,HPX::Tuple, spec_ind::Vector{Int})
        X=zeros(length(gas.X))
        X[spec_ind[1:end-1]].=HPX[end]
        X[spec_ind[end]] = 1-sum(X)
        setHPX(gas,(HPX[1],HPX[2],X))
        return nothing
    end


    function setUVX(gas::PyObject, UVX::Tuple)
        gas.UVX=UVX
    end

    function int_nrg(gas::PyObject)
        return gas.int_energy_mass
    end

    function enthalpy(gas::PyObject)
        return gas.enthalpy_mass
    end

    function cv(gas::PyObject)
        return gas.cv
    end

    function cp(gas::PyObject)
        return gas.cp
    end

    function rho(gas::PyObject)
        return gas.density
    end
    
    function P(gas::PyObject)
        return gas.P
    end

    function T(gas::PyObject)
        return gas.T
    end

    function MW_spec(gas::PyObject)
        return gas.molecular_weights
    end

    function MW_mix(gas::PyObject)
        return gas.mean_molecular_weight
    end

    function spec_inds(gas::PyObject, specs::Vector{String})
        return stack([gas.species_index(s) for s in specs])
    end

    function Y(gas::PyObject)
        return gas.Y
    end

    function X(gas::PyObject)
        return gas.X
    end

    function u0(gas::PyObject,specs::Vector{String})
        u=Vector{Float64}(undef,length(specs))
        for i in eachindex(specs)
            s=specs[i]
            gas.TPX=(298.15,101.325e3,"$s:1.0")
            u[i]=int_nrg(gas)
        end
        return u
    end

    function h0(gas::PyObject,specs::Vector{String})
        u=Vector{Float64}(undef,length(specs))
        for i in eachindex(specs)
            s=specs[i]
            gas.TPX=(298.15,101.325e3,"$s:1.0")
            u[i]=enthalpy(gas)
        end
        return u
    end

end