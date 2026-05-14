function LHV_calc(Hf, comp)

end

function Hf_gas(gas, spec::String; mass_basis=true)
    gas.TPX=(298.15, 101.325e3, @sprintf("%s:1.0", spec))
    mass_basis ? hf=gas.enthalpy_mass : hf=gas.enthalpy_mole
    return hf
end

function Hf_gas(gas, spec::Vector{String}; mass_basis=true)
    Hf=stack([Hf_gas(gas,s, mass_basis=mass_basis) for s in spec]) # in J/kmol
    return Hf
end

function get_LHV(gas, specs)
    Hf=Hf_gas(gas, specs, mass_basis=false)/1e6 # in kJ/mol
    LHV=get_LHV(gas,specs,Hf)
    return LHV
end

function get_LHV(gas, specs, Hf)
    # assume every carbon goes to CO2
    # every hydrogen goes to water
    # ignore everything else
    # all species gaseous
    Hf_CO2=Hf_gas(gas,"CO2", mass_basis=false)/1e6
    Hf_H2O=Hf_gas(gas,"H2O", mass_basis=false)/1e6
    LHV=Vector{Float64}(undef, length(specs))
    for i in eachindex(Hf)
        LHV[i] = Hf_CO2 * gas.n_atoms(specs[i],"C") + 
            Hf_H2O * gas.n_atoms(specs[i],"H")/2 - 
            Hf[i]
    end
    return -LHV
end
