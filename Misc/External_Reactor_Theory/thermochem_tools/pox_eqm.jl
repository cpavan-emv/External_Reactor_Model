function set_oxidizer(gas,Fspec::Vector{String}, Fcomp::Vector{Float64},ϕ, xO2)
    cnt=[0.0,0.0,0.0] # count of C, H and O in fuel
    atom_names=["C","H","O"]
    Fcomp /=sum(Fcomp)
    for i in eachindex(Fspec)
        atoms=collect(keys(gas.species(Fspec[i]).composition))
        for j in 1:3
            if atom_names[j] in atoms
                cnt[j] += gas.species(Fspec[i]).composition[atom_names[j]] * Fcomp[i]
            end
        end
    end
    f_o_stoich=1/(cnt[1]+cnt[2]/4-cnt[3]/2)
    f_o_actual=f_o_stoich*ϕ
    O2_actual=1/f_o_actual
    N2_actual=(1-xO2)/xO2*O2_actual
    if any(Fspec .== "O2")
        Fcomp[Fspec .== "O2"] .+= O2_actual
    else
        Fspec=vcat(Fspec,"O2")
        Fcomp=vcat(Fcomp,O2_actual)
    end 
    if any(Fspec .== "N2")
        Fcomp[Fspec .== "N2"] .+= N2_actual
    else
        Fspec=vcat(Fspec,"N2")
        Fcomp=vcat(Fcomp,N2_actual)
    end 
    return comp_string(Fcomp,Fspec), N2_actual+O2_actual
end

function find_pox_eqm(gas, ϕ, xO2, fuel_specs, fuel_comp; 
    Teq=1600.0, Peq=100e3, specs_out=["CO","CO2","H2","H2O","N2","Other"])
    #=
    Calculates POX equilibrium 
    returns the mole fraction of "specs out" and the molar flowrate out
    per mole of fuel in
    =#
    comp, ndot_ox=set_oxidizer(gas, fuel_specs, fuel_comp, ϕ, xO2)
    gas.TPX = (Teq, Peq, comp)
    MW0=gas.mean_molecular_weight
    gas.equilibrate("TP")
    MW1=gas.mean_molecular_weight
    X=[gas.X[gas.species_index(s)+1] for s in specs_out[1:end-1]]
    nf_no=1/(1+ndot_ox) * MW1/MW0 # ratio of fuel flow to outlet flow
    return vcat(X, 1-sum(X)), 1/nf_no
end

