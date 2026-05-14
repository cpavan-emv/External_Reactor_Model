# this version moves all DE tools to separate 

# first get the environment set up
topdir=abspath((@__DIR__)*"/../")
using Pkg
Pkg.activate(topdir*"/src/.")
using LinearAlgebra

include(topdir*"src/fluid_properties/fluid_props.jl")
include(topdir*"src/DE_Construction/DE_Const.jl")
include(topdir*"src/evaluation/conservation_checks.jl")
DE_construct.fluid_props=fluid_props

# initialize the gas objects and user-set parameters
gas1=fluid_props.initialize_ideal_gas("gri30.yaml")
gas2=fluid_props.initialize_ideal_gas("gri30.yaml")
fluid_props.kingas=fluid_props.initialize_ideal_gas("gri30.yaml")
rhocat=1e3 # kg/m^3

specs=["CO", "CO2", "H2", "CH3OH", "H2O", "N2"]
DE_construct.set_gas_constants(gas1,specs)

# Define the compression function

function V_dV(t, CR, tcomp)
    # assume that V follows a cosine curve
    V=(CR-1)/2*(cos(pi*t/tcomp)+1)+1e-3
    dV=(CR-1)/2*(-1)*pi/tcomp*sin(pi*t/tcomp)
    if t>tcomp
        V=1e-3
        dV=0.0
    end
    return V, dV

end

# setup the initial conditions and the parameters
CR=8
P0=1
tcomp=10e-3
Vfunc(t)=V_dV(t, CR, tcomp)
gas_props0=(300.0, 100e3, "H2:0.01, CO2:0.01, N2:0.98")
#gas_props0=(300.0, 100e3, "N2:1.0")
fluid_props.setTPX(gas1, gas_props0)
fluid_props.setTPX(gas2, gas_props0)
param=(gas1, gas2, Vfunc, rhocat, 250+273.15)

spec_ind=DE_construct.spec_ind
u0=vcat([P0*100e3, 300.0, 300.0], 
    fluid_props.X(gas1)[spec_ind[1:end-1]], 
    fluid_props.X(gas2)[spec_ind[1:end-1]])
du=copy(u0)

using OrdinaryDiffEq
tmax=2*tcomp;
tp=range(0.0,tmax,201) # times for saving solution

# full implicit solution
prob=ODEProblem(DE_construct.f!,u0,(0.0,tmax))
sol = solve(prob,RadauIIA5(autodiff=false),p=param, reltol=1e-8, abstol=1e-10, saveat=tp)#
y=permutedims(stack(sol.u))

V=permutedims(stack(Vfunc.(tp)))[:,1]
dV=permutedims(stack(Vfunc.(tp)))[:,2]

######################################
# load RCM video
using VideoIO
file=topdir*"data/videos/RCM_ASM_10mmps_10fps_1mmpf.mp4"

# each frame in this video moves the RCM piston 1mm
# roughly 200mm corresponds to a CR of 8
# stroke = 200mm, CR = 8 => CR=(Vf+Vs)/Vf=1+Vs/Vf=>Vf=Vs/(CR-1)
Vf=200/(CR-1) 
# convert these to a corresponding volume ratio
# Vstar = Vs/Vf
Nf=VideoIO.get_number_frames(file)
frame_V=collect((200 .- (0:Nf))/Vf)
fnum=[argmin(abs.(vtmp .- frame_V)) for vtmp in V]


########################################

using Plots, Measures, LaTeXStrings 
begin 
global ind_read=1
global ind_read_last=1
global cnt_read=0
vid=VideoIO.openvideo(file)
img=read(vid)
anim =@animate for i in eachindex(tp)
    global cnt_read
    if i==1
        global ind_read_last=1
    end
    p=[Plots.plot() for _ in 1:2]
    while fnum[i]>ind_read
        global ind_read+=1
        EOF=skipframe(vid, throwEOF=false)
        cnt_read+=1
    end
    
    p[2]=plot(tp[1:i]*1e3,y[1:i,1]/100e3, 
        xlabel="Time (ms)", ylabel="Pressure (bar)",ylim=(0,20), xlim=(0, tp[end]*1e3),
        label="Pressure", leftmargin=10mm,bottommargin=10mm)
    plot!(p[2], [NaN], [NaN], color=:black,label="Volume ratio")
    plot!(twinx(p[2]),tp[1:i]*1e3, (V[1] .+ 1) ./ (V[1:i] .+ 1), ylabel=L"V_0/V",
        ylim=(0,10), xlim=(0, tp[end]*1e3), color=:black, legend=:none, rightmargin=10mm)
    
    if ind_read!=ind_read_last && ind_read<Nf
        read!(vid, img);
        global ind_read+=1
    end
    p[1]=plot(img[:,900:1760], xaxis=false, yaxis=false)
    fig=plot(p..., size=(1080,720))
    global ind_read_last = copy(ind_read)
end
VideoIO.close(vid)
gif(anim, "RCM_demo.gif", fps=10)
end



savefig(fig,topdir*"output/figures/heat_loss_validation.png")


DF=DataFrame(:t=>tp, :P=>y[:,1])
CSV.write("heat_loss_validation_data_segmented.csv", DF)

# #####################################################################
# # Conservation checks
# yeval=y
# rho_uA=get_rhoU(gas1,yeval,1)
# rho_uB=get_rhoU(gas2,yeval,2)
# M=Mtot(rho_uA[:,1], rho_uB[:,1], V)
# plot(tp,M)

# U=Utot(rho_uA,rho_uB,V)
# Ein=Ecomp(yeval[:,1],V)
# plot(tp,U-Ein,legend=:bottomright)