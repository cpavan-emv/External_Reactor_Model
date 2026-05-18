# first get the environment set up
topdir = abspath((@__DIR__)*"/../../")  # → project root
using Pkg
Pkg.activate(topdir*"src/.")
using LinearAlgebra

include(topdir*"src/DE_Construction/DE_Model.jl")
include(topdir*"src/Engine_Functions.jl")
include(topdir*"src/conservation_checks.jl")
const ct = DE_Model.ct
const mech_file=abspath(topdir*"data/mechanisms/MeOH_Synth.yaml")

###########################################################
# Allocate cases
beta_array=collect(0.5:0.25:5)
CA_open_array=(collect(10:2:30) .* ones(1,length(beta_array)))'
beta_array=beta_array .* ones(size(CA_open_array,2))'
cases=Matrix{Tuple{Float64,Int}}(undef, size(CA_open_array))
for i in axes(beta_array,1), j in axes(beta_array,2)
    cases[i,j]=(beta_array[i,j], CA_open_array[i,j])
end
cases=Vector(cases[:])

Nproc=parse(Int, ARGS[1])
myind=parse(Int,ARGS[2]) # indexes from 1 to Nproc

Ncase_each=Int(floor(length(cases)/Nproc))
myrng=((myind-1)*Ncase_each+1):min(myind*Ncase_each,length(cases))

mycases=cases[myrng]
# also assign extra cases
rem=length(cases)-Nproc*Ncase_each
if myind<=rem
    push!(mycases, cases[Nproc*Ncase_each+myind])
end


###########################################################
gases = [DE_Model.initialize_ideal_gas(mech_file) for _ in 1:3]
rhocat = 1e3  # kg/m^3 catalyst effective density

# Volume definitions
# V = displaced volume / clearance volume; V→0 at TDC
Vdisp_dimensional  = 3e-4   # piston displacement (reference dimension)
CR0    = 18  # compression ratio with valve closed (CR = (Vdisp+Vclear)/Vclear)
Vclear = Vdisp_dimensional / (CR0 - 1)

# add in engine dimensions
R=3.5
RPM=600
ω=RPM/60*2*π
Vclear=Vdisp_dimensional*1/(CR0-1)

function V_dV(t)
    V, dV = V_dVfunc(t, ω, CR0,R) # this calculates Vtotal/Vclearance
    V *= Vclear/Vdisp_dimensional # this is Vtotal/Vdisp
    dV *= Vclear/Vdisp_dimensional
    return V, dV
end

# Initial conditions and parameters
Pin    = 1   # bar (intake manifold)
Pout   = 1  # bar (exhaust / reactor side)
Preact = 10  # bar (external reactor starting pressure)

Vfunc=V_dV
gas_props0 = (250+273.15, 100e3, "H2:0.75, CO2:0.25, CO:0.0")
gas_props1 = (300.0, 100e3, "N2:1.0")
foreach(g -> ct.setTPX(g, gas_props1), gases[1:2])

Nspec = gases[1].gas.Nspec

u0_ext = [[Preact*100e3, 250+273.15]; gases[1].X[1:end-1]]
u0_cyl = [[Pin*100e3, 25+273.15]; gases[2].X[1:end-1]]
u0     = [u0_ext; u0_cyl]

gas_tmp = gases[1]
ct.setTPX(gas_tmp, gas_props0)

TPX_intake  = (25+273.15, Pin*100e3,  gas_tmp.X)
TPX_exhaust = (25+273.15, Pout*100e3, gas_tmp.X)

# Mode objects — created once, reused every cycle
mode_ex      = DE_Model.IntakeExhaust(TPX_exhaust)
mode_in      = DE_Model.IntakeExhaust(TPX_intake)
mode_coupled = DE_Model.Coupled()

trev  = 1/(RPM/60)
tcomp=trev/2
tp  = range(0.0, trev, 1001)
tex = tp[1:floor(Int, length(tp)/2)+1]
tin = tp[floor(Int, length(tp)/2)+1:end]

##################
Ncycle = 10
tcycle = trev*2


import BSON, DataFrames
using Printf
###################################
for i in eachindex(mycases)
    @printf("Beginning case %u on worker %u\n", i, myind)
    c=mycases[i]

CA_open=c[2] # degrees
beta    = c[1] # external volume as a multiple of displacement volume
params = DE_Model.ReactorParams(gases, Vfunc, beta, Vdisp_dimensional, rhocat,
                                [250.0, 85.0, 85.0] .+ 273.15,
                                [50, 500, 500]*1e-3,
                                1e-1*[2.0, 1, 6.0, 10.0])  # Cv [Cv_BA_comp, Cv_BA_exp, Cv_IE_exh, Cv_IE_int]

topen = CA_open * π/180/ω
u = u0'
t = 0.0

for N in 1:Ncycle
    tstart = t[end]

    # BDC → exhaust stroke
    tloc, y = DE_Model.evolve(mode_ex, u[end, :], params, tex)
    t = [t; tloc .+ tstart]
    u = [u; y]

    # Intake stroke
    tloc, y = DE_Model.evolve(mode_in, u[end, :], params, tin)
    t = [t; tloc .+ tstart]
    u = [u; y]
    tstart = t[end]

    # Compression — decoupled until pressures equalise
    cond(u, t, int) = DE_Model.Vars_Eq_condition(u, t, int, [1, Nspec+2])
    tloc, y = DE_Model.evolve(DE_Model.Decoupled(), u[end, :], params, tp, cond)
    t = [t; tloc .+ tstart]
    u = [u; y]

    # Coupled phase (chambers connected) through TDC + topen
    tp2 = [tloc[end]; tp[tp .> tloc[end]]]
    tp2 = [tp2[tp2 .< tcomp+topen]; tcomp+topen]
    tloc, y = DE_Model.evolve(mode_coupled, u[end, :], params, tp2)
    t = [t; tloc .+ tstart]
    u = [u; y]

    # Expansion — decoupled again
    tp3 = [tp2[end]; tp[tp .> tp2[end]]]
    tloc, y = DE_Model.evolve(DE_Model.Decoupled(), u[end, :], params, tp3)
    t = [t; tloc .+ tstart]
    u = [u; y]

    println("Cycle $N Complete")
end

file=@sprintf("output/data/topen=%udeg_Vc=%.1fVd.bson",Int(CA_open), beta)
BSON.bson(file,Dict(:t=>t, :u=>u, :case=>c))
end