const Ru=8.314 # J/mol K

function G(T,P,spec)
    kingas.TPX=T,P, @sprintf("%s:1", spec)
    return kingas.enthalpy_mole - kingas.entropy_mole * T
end

# equilibrium constants
function K_CO2hyd(T,P)
    dG=G(T,P,"CH3OH")+G(T,P,"H2O")-
        G(T,P,"CO2")-3*G(T,P,"H2")
    return exp(-dG/(Ru*1e3*T))*(100e3/P)^2 # bar^-2
end
function K_COhyd(T,P)
    dG=G(T,P,"CH3OH")-
        G(T,P,"CO")-2*G(T,P,"H2")
    return exp(-dG/(Ru*1e3*T))*(100e3/P)^2 # bar^-2
end
function K_RWGS(T,P)
    dG=-G(T,P,"CO2")-G(T,P,"H2") + 
        G(T,P,"CO")+G(T,P,"H2O")
    return exp(-dG/(Ru*1e3*T))
end

# rates from Bisotti
# adsorption constants (Table 7)
# or-GR
bCO(T)=2.16e-5*exp(46800/(Ru*T)) # bar ^-1
bCO2(T)=7.05e-7*exp(61700/(Ru*T)) # bar ^-1
bH2O_H2(T)=6.37e-9*exp(84000/(Ru*T)) # bar ^-1/2
# ref-GR
# bCO(T) = 1.540e-3*exp(14936/(Ru*T)) # bar ^-1
# bCO2(T) = 8.206e-9*exp(76594/(Ru*T)) # bar ^-1
# bH2O_H2(T) = 3.818e-9*exp(97350/(Ru*T)) # bar ^-1/2

# kinetic rates (Table 8)
# or-GR
k_COhyd(T) = 4.89e7*exp(-113000/(Ru*T)) # mol/s/bar/kgcat
k_RWGS(T) = 9.64e11*exp(-152900/(Ru*T)) # mol/s/bar^(1/2)/kgcat
k_CO2hyd(T) = 1.09e5*exp(-87500/(Ru*T)) # mol/s/bar/kgcat
# ref-GR
# k_COhyd(T) = 2.240e7*exp(-106729/(Ru*T)) # mol/s/bar/kgcat
# k_RWGS(T) = 4.241e13*exp(-149856/(Ru*T)) # mol/s/bar^(1/2)/kgcat
# k_CO2hyd(T) = 9.205e1*exp(-45889/(Ru*T)) # mol/s/bar/kgcat


# reaction rates
# x will be mole fractions ordered as [CO, CO2, H2, CH3OH, H2O]
# units should be mol/s/kgcat when P is in bar
denom(T, P, x)=(1 .+ bCO(T)*P*x[1]+bCO2(T)*P*x[2]).*
    ((P*x[3]).^(1/2)+(bH2O_H2(T))*P*x[5])

r1(T,P,x) = k_COhyd(T)*bCO(T)*(
    (P*x[1].*(P*x[3]).^(3/2)-P*x[4]./((P*x[3]).^(1/2)*K_COhyd(T,100e3)))./
    denom(T,P,x))
r2(T,P,x) = k_RWGS(T)*bCO2(T)*(
        (P*x[2].*P.*x[3]-(P*x[1].*P.*x[5])./(K_RWGS(T,100e3)))./
        denom(T,P,x))
r3(T,P,x) = k_CO2hyd(T)*bCO2(T)*(
    (P*x[2].*(P*x[3]).^(3/2)-P*x[4]*P.*x[5]./((P*x[3]).^(3/2)*K_CO2hyd(T,100e3)))./
    denom(T,P,x))

    