using LinearAlgebra, Polynomials, DataFrames, CSV
# calculates real gas properties using the SRK equation of state
# See Bennekom et al dx.doi.org/10.1021/ie3017362 | Ind. Eng. Chem. Res. 2012, 51, 12233−12243
# See also ref 20 therin (Ind. Eng. Chem. Process Des. Dev., Vol. 18, No. 2, 1979)
# for brief summary of method

# Binary interaction coefficients
# taken from Ind. Eng. Chem. Res. 2022, 61, 2206−2226
# https://doi.org/10.1021/acs.iecr.1c04476
# order is CO, CO2, H2, H2O, CH4, N2, O2, CH3OH

# const BIC=[0    0.1164  -7e-4   -0.5594 2.04e-2 1.30e-2  0   0;
#     0.1164  0   0.1164  -0.12155    9.56e-2  -1.71e-2  9.75e-2  1.7e-2;
#     -7e-4   0.1164  0   -0.7544 1e-4    -1e-3   0   0;
#     -0.5594 -0.12155    -0.7544 0   0.5 -0.69648    0   -9e-2;
#     2.04e-2 9.56e-2 1e-4    0.5 0   3.12e-2 0   -3.5e-2;
#     1.3e-2  -1.71e-2    -1e-3   -0.69648    3.12e-2 1   -1.4e-2 -0.2141;
#     0   9.75e-2 0   0   0   -1.4e-2 0 0;
#     0   1.7e-2  0   -9e-2   -3.5e-2 -0.2141 0   0]

# taken from Bennekom supplementary material. O2, N2 parametes from above
#           CO      CO2         H2          H2O         CH4         N2          O2          CH3OH
const BIC=[ 0       0.1164      -7e-4       -0.474      2.04e-2     1.30e-2     0           -0.37;
            0.1164  0           0.1164      0.3         0.0956      -1.71e-2    9.75e-2     0.1;
            -7e-4   0.1164      0           -0.745      1e-2        -1e-3       0           -0.125;
            -0.474  0.3         -0.745      0           0.014       -0.69648    0           -0.075;
            2.04e-2 9.56e-2     0.001       0.014       0           3.12e-2     0           0.046;
            1.3e-2  -1.71e-2    -1e-3       -0.69648    3.12e-2     0           -1.4e-2     -0.2141;
            0       9.75e-2     0           0           0           -1.4e-2     0           0;
            -0.37   0.1         -0.125      -0.075      0.046       -0.2141     0           0]


# from Aspen
const Pc=[34.99; 73.83; 13.13; 220.64; 45.99; 34.0; 50.43; 80.84]*100e3
const Tc=[134.45; 304.2; 33.18; 647.1; 190.6; 126.2; 154.58; 513.0]
const ω=[0.039465, 0.22551, -0.22014, 0.34417, 0.011433, 0.036816, 0.021324, 0.56197]

begin
    specs=["CO", "CO2", "H2","H2O","CH4","N2","O2","CH3OH"];
    const SRK_spec=Dict(specs[i] => i for i in eachindex(specs))
    specind(spec::String) = SRK_spec[spec]
end

# these are all valid for pure fluids
a(spec::Int)=0.42748*Ru^2*Tc[spec]^2/Pc[spec] # There is a typo in Bennekom paper - should be Pc^1 not Pc^2
a(spec::String)=a(SRK_spec[spec])

b(spec::Int)=0.08664*Ru*Tc[spec]/Pc[spec]
b(spec::String)=b(SRK_spec[spec])

# I am using the standard equation (wikipedia) rather than the Bennekom one
# wikipedia references the Graboski and Dauber paper (ref 20 in Bennekom)
# Note special formulation for H2 (makes minimal difference in my system)
α(spec::Int,T::Float64)=spec==3 ? 1.202*exp(-0.30288*T/Tc[spec]) : (1+(0.48508+1.55171*ω[spec]-0.15613*ω[spec]^2)*(1-sqrt(T/Tc[spec])))^2
α(spec::String, T::Float64)=α(SRK_spec[spec],T)


function Z(p,T, spec::String)
    # pure species compressibility factor found from cubic EoS (SRK)
    A=a(spec)*α(spec,T)*p / (Ru* T)^2
    B=b(spec)*p/(Ru*T)
    terms=[-A*B, # x0
        A-B-B^2,# x1
        -1, #x2
        1] # x3
    return maximum(real.(roots(Polynomial(terms))))
end

function fugacity_coeff(p,T,spec::String)
    # pure species fugacity coeff
    Zcalc=Z(p,T,spec)
    Acalc=a(spec)*α(spec,T)*p / (Ru* T)^2
    Bcalc=b(spec)*p/(Ru*T)
    return exp(Zcalc-1.0-log(Zcalc - Bcalc)-Acalc/Bcalc*log(1+Bcalc/Zcalc))
end

# these are for mixtures
function abmix(T::Float64)
    # first construct full BIC matrix for a
    amat=diagm(map(i -> a(i)*α(i,T),eachindex(Pc)))
    for i in range(1, size(amat,2)-1), j in range(i+1, size(amat,2))
        amat[i,j]=sqrt(amat[i,i]*amat[j,j])*(1-BIC[i,j])
        amat[j,i]=amat[i,j]
    end
    bmat=diagm(map(i -> b(i),eachindex(Pc)))
    for i in range(1, size(amat,2)-1), j in range(i+1, size(amat,2))
        bmat[i,j]=(bmat[i,i]+bmat[j,j])/2
        bmat[j,i]=bmat[i,j]
    end
    return amat, bmat
end

function construct_xvec(x::Vector{Float64}, specs::Vector{String})
    x1=zeros(Float64, length(Pc))
    x/=sum(x)
    for i in eachindex(specs)
        x1[SRK_spec[specs[i]]]=x[i]
    end
    return x1
end

function Z(p,T, x::Vector{Float64})
    # pure species compressibility factor found from cubic EoS (SRK)
    am,bm=abmix(T)
    amix=sum(am .* (x' .*x))
    bmix=sum(bm .* (x' .*x))
    Acalc=amix*p/(Ru^2 * T^2)
    Bcalc=bmix*p/(Ru*T)
    terms=[-Acalc*Bcalc, # x0
        Acalc-Bcalc-Bcalc^2,# x1
        -1, #x2
        1] # x3
    return maximum(real.(roots(Polynomial(terms))))
end

function fugacity_coeff(p,T,x::Vector{Float64})
    # partial molar fugacities
    Zmix=Z(p,T,x)
    am,bm=abmix(T)
    amix=sum(am .* (x*x'))
    bmix=sum(bm .* (x*x'))
    #bmix=sum([x[i]*b(i) for i in eachindex(x)])

    Amix=amix*p/(Ru^2 * T^2)
    Bmix=bmix*p/(Ru*T)

    ϕ=[exp(b(i)*(Zmix-1)/bmix - 
        log(Zmix-Bmix)-
        Amix/Bmix*(
            2*sum(x .* am[i,:])/amix-b(i)/bmix
        )*log(1 .+ Bmix/Zmix)) for i in eachindex(x)]

    return ϕ
end

function fugacity_coeff_mix(p,T,x)
    # fugacity of the one-fluid mixture
    am,bm=abmix(T)
    amix=sum(am .* (x' .*x))
    bmix=sum(bm .* (x*x'))
    #bmix=sum([x[i]*b(i) for i in eachindex(x)])
    Zmix=Z(p,T,x)
    Amix=amix*p/(Ru^2 * T^2)
    Bmix=bmix*p/(Ru*T)
    return exp(Zmix-1.0-log(Zmix - Bmix)-Amix/Bmix*log(1+Bmix/Zmix))
end

# load in the WGS and CO hydrogenation equilibrium constants calculated by Aspen
function load_Keq(dat_fil)
    df=CSV.read(dat_fil, DataFrame)
    thermochem_tools.TK=df.T
    thermochem_tools.KWGS=df.K_WGS
    thermochem_tools.KCOh=df.K_Cohyd
    return nothing
end