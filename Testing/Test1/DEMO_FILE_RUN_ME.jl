# first get the environment set up
topdir = abspath((@__DIR__)*"/../../")  # → project root
using Pkg
Pkg.activate(topdir*"src/.")
using LinearAlgebra

include(topdir*"src/DE_Construction/DE_Model.jl")
include(topdir*"src/conservation_checks.jl")
const ct = DE_Model.ct
const mech_file=abspath(topdir*"data/mechanisms/MeOH_Synth.yaml")

gases = [DE_Model.initialize_ideal_gas(mech_file) for _ in 1:3]
rhocat = 1e3  # kg/m^3 catalyst effective density


# Volume definitions
# V = displaced volume / clearance volume; V→0 at TDC
Vdisp_dimensional  = 3e-4   # piston displacement (reference dimension)
beta    = 2 # external volume as a multiple of displacement volume
CR0    = 18  # compression ratio with valve closed (CR = (Vdisp+Vclear)/Vclear)
Vclear = Vdisp_dimensional / (CR0 - 1)

function V_dV(t, CR, param)
    tcomp  = param[1]
    tdelay = param[2]
    V  =  1/2 * (cos(pi*(t-tdelay)/tcomp)+1) + 1/CR
    dV = -1/2 * pi/tcomp * sin(pi*(t-tdelay)/tcomp)
    return V, dV
end

# Initial conditions and parameters
Pin    = 4   # bar (intake manifold)
Pout   = 10  # bar (exhaust / reactor side)
Preact = 10  # bar (external reactor pressure)

tcomp  = 50e-3
tdelay = 0
Vfunc(t) = V_dV(t, CR0, (tcomp, tdelay))

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

params = DE_Model.ReactorParams(gases, Vfunc, beta, Vdisp_dimensional, rhocat,
                                [250.0, 85.0, 85.0] .+ 273.15,
                                [50, 500, 500]*1e-3,
                                1e-1*[2.0, 1, 6.0, 10.0])  # Cv [Cv_BA_comp, Cv_BA_exp, Cv_IE_exh, Cv_IE_int]

# Mode objects — created once, reused every cycle
mode_ex      = DE_Model.IntakeExhaust(TPX_exhaust)
mode_in      = DE_Model.IntakeExhaust(TPX_intake)
mode_coupled = DE_Model.Coupled(Float64(beta))

tmax  = tcomp * 2
topen = tcomp / 15
tp  = range(0.0, tmax, 1001)
tex = tp[1:floor(Int, length(tp)/2)+1]
tin = tp[floor(Int, length(tp)/2)+1:end]

##################
Ncycle = 100
tcycle = tcomp * 2

u = u0'
t = 0.0

for N in 1:Ncycle
    global t, u
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

yA = u[:, 1:1+Nspec]
yB = u[:, Nspec+2:2*(1+Nspec)]

specs = ct.get_speciesName_all(gases[1].gas)

using Plots, Measures
begin
    plt = [Plots.plot() for _ in 1:4]
    plot!(plt[1], t, [yA[:,2] yB[:,2]] .- 273,
        xlabel="Time (s)", ylabel="Temp (degC)",
        label=["Ext. Vol" "Cyl Vol"])
    plot!(plt[2], t, [yA[:,1] yB[:,1]]/100e3,
        xlabel="Time (s)", ylabel="Pressure (bar)",
        label=["Ext. Vol" "Cyl Vol"])
    XA = yA[:, 3:end]
    XA = [XA  1 .- sum(XA, dims=2)]
    XB = yB[:, 3:end]
    XB = [XB  1 .- sum(XB, dims=2)]
    plot!(plt[3], t, XA*100, xlabel="Time (s)", ylabel="Mole Fraction (ext. vol)",
        ylim=(0,100), label=permutedims(specs))
    plot!(plt[4], t, XB*100, xlabel="Time (s)", ylabel="Mole Fraction (Cyl vol)",
        ylim=(0,100), label=permutedims(specs))
    for p in plt
        plot!(p, xlim=(4*(Ncycle-2)*tcomp, 4*Ncycle*tcomp))
        #plot!(p, xlim=(0*tcomp, 4*Ncycle*tcomp))
    end
    plot(plt..., size=(1080, 720), margin=10mm)
end
