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

gases = [DE_Model.initialize_ideal_gas(mech_file) for _ in 1:3]
rhocat = 1e3  # kg/m^3 catalyst effective density

# Volume definitions
# V = displaced volume / clearance volume; V→0 at TDC
Vdisp_dimensional  = 3e-4   # piston displacement (reference dimension)
beta    = 2 # external volume as a multiple of displacement volume
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

params = DE_Model.ReactorParams(gases, Vfunc, beta, Vdisp_dimensional, rhocat,
                                [250.0, 85.0, 85.0] .+ 273.15,
                                [50, 500, 500]*1e-3,
                                1e-1*[2.0, 1, 6.0, 10.0])  # Cv [Cv_BA_comp, Cv_BA_exp, Cv_IE_exh, Cv_IE_int]

# Mode objects — created once, reused every cycle
mode_ex      = DE_Model.IntakeExhaust(TPX_exhaust)
mode_in      = DE_Model.IntakeExhaust(TPX_intake)
mode_coupled = DE_Model.Coupled()

trev  = 1/(RPM/60)
tcomp=trev/2
CA_open=18 # degrees
topen = CA_open * π/180/ω
tp  = range(0.0, trev, 1001)
tex = tp[1:floor(Int, length(tp)/2)+1]
tin = tp[floor(Int, length(tp)/2)+1:end]

##################
data_folder="/mnt/c/Users/cpava/Downloads/ext_reactor_pressure/data1/"
files=readdir(data_folder)
import BSON, DataFrames

P_P0_all=Matrix{Float64}(undef, Ncycle, length(files))
dM_frac_all=similar(P_P0_all)
dM=similar(P_P0_all)
cases=Vector{Tuple{Float64,Int}}(undef, length(files))

for i in eachindex(files)

file=joinpath(data_folder,files[i])
result = BSON.load(file)
cases[i]=result[:case]


u=result[:u];
t=result[:t];
Ncycle=250


yA = u[:, 1:1+Nspec]
yB = u[:, Nspec+2:2*(1+Nspec)]

specs = ct.get_speciesName_all(gases[1].gas)

# using Plots, Measures
# begin
#     plt = [Plots.plot() for _ in 1:4]
#     plot!(plt[1], t, [yA[:,2] yB[:,2]] .- 273,
#         xlabel="Time (s)", ylabel="Temp (degC)",
#         label=["Ext. Vol" "Cyl Vol"])
#     plot!(plt[2], t, [yA[:,1] yB[:,1]]/100e3,
#         xlabel="Time (s)", ylabel="Pressure (bar)",
#         label=["Ext. Vol" "Cyl Vol"])
#     XA = yA[:, 3:end]
#     XA = [XA  1 .- sum(XA, dims=2)]
#     XB = yB[:, 3:end]
#     XB = [XB  1 .- sum(XB, dims=2)]
#     plot!(plt[3], t, XA*100, xlabel="Time (s)", ylabel="Mole Fraction (ext. vol)",
#         ylim=(0,100), label=permutedims(specs))
#     plot!(plt[4], t, XB*100, xlabel="Time (s)", ylabel="Mole Fraction (Cyl vol)",
#         ylim=(0,100), label=permutedims(specs))
#     for p in plt
#         #plot!(p, xlim=(4*(Ncycle-2)*tcomp, 4*Ncycle*tcomp))
#         plot!(p, xlim=(0*tcomp, 4*Ncycle*tcomp))
#     end
#     plot(plt..., size=(1080, 720), margin=10mm)
# end

rhoU_A = get_rhoU(gases[1], yA, 1)
rhoU_B = get_rhoU(gases[2], yB, 1)
V=stack(Vfunc.(t))[1,:]

dM_in=zeros(Float64,Ncycle)
dM_out=zeros(Float64,Ncycle)
dM_frac=similar(dM_in)
P_P0=similar(dM_in)
using Printf

ind=length(t)
for N in Ncycle:-1:1
    # time cycle ends
    tfinal=t[end]-(Ncycle-N)*(2*trev)
    while t[ind]>tfinal
        ind-=1
    end
    Mcyl_o=rhoU_B[ind]*V[ind] # mass in cylinder at time outlet opens
    Mext_o=rhoU_A[ind]*beta

    # find ind for TDC compression
    while t[ind]>tfinal-tcomp
        ind-=1
    end
    Mcyl_f=rhoU_B[ind]*V[ind] # mass in cylinder at TDC
    Mext_f=rhoU_A[ind]*beta

    # find a time before valve to external reactor opens
    while t[ind]>tfinal-tcomp-topen*2
        ind-=1
    end
    # grab external volume pressure here
    P_P0[N]=yA[ind,1]/(Pin*100e3)
    Mcyl_i=rhoU_B[ind]*V[ind] # mass in cylinder before valve opents
    Mext_i=rhoU_A[ind]*beta

    # conservation check - these should be equal
    dM_cyl=Mcyl_i-Mcyl_f # mass moving to external volume (kg/cycle/Vdisp)
    dM_ext=Mext_f-Mext_i

    #@printf("Mass Imbalance on cycle %u: %.2e\n",N,dM_cyl-dM_ext)
 
    dM_in[N]=dM_cyl
    dM_out[N]=Mcyl_o-Mcyl_f
    dM_frac[N]=(Mcyl_i-Mcyl_f)/Mcyl_i # fraction of the mass in cylinder moved to external volume
end
    
P_P0_all[:,i] .= P_P0
dM_frac_all[:,i] .= dM_frac
dM[:,i] .= dM_in

end

case_name=[@sprintf("beta=%.1f, CA_open=%u", c[1], c[2]) for c in cases]

using Plots, Measures
plot(P_P0_all, label=permutedims(case_name), 
    xlabel="Cycle Number", ylabel="P/P_0",
    palette=:tab20, size=(1080,720),
    margin=10mm)

plot(dM_frac_all, label=permutedims(case_name), 
    xlabel="Cycle Number", ylabel="Pumping Efficiency",
    palette=:tab20, size=(1080,720),
    margin=10mm)

plot(dM*RPM/120*Vdisp_dimensional*1e3, label=permutedims(case_name), 
    xlabel="Cycle Number", ylabel="mass flow (g/s)",
    palette=:tab20, size=(1080,720),
    margin=10mm)