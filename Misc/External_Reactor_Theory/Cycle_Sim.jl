topdir=abspath((@__DIR__)*"./../../")
using Pkg
Pkg.activate(topdir*"Kinetic_Model/.")

#################################
CR=20
RPM=1200
ω=RPM*2*pi/60
γ=1.4
L=82e-3
ℓ=137e-3
τ=2*π/ω # cycle time (1 rev, not 1 engine cycle)
Sp0=L/(τ/2) # average piston speed
R=ℓ/(L/2)
rc=CR
Vd=1.0e-3/3
Vc=Vd/(rc-1)
A=Vd/L


function Vfun0(θ)
    #if abs(θ) .> π
    #    return 1
    #else
    return (1 .+ 0.5*(rc-1) .* (R .+ 1.0 .- cos.(θ) .- ((R^2 .- (sin.(θ)).^2).^(1/2))))*Vc/(Vc+Vd)
    #end
end

CA=-360:0.5:360
CA_rad=CA*π/180



###########################
# Baseline

V0=Vfun0(CA_rad)
P0=V0.^(-γ)
P0[abs.(CA_rad) .>= π] .= 1.0
using Plots, LaTeXStrings
begin
    plt=[Plots.plot() for _ in 1:2]
    plot!(plt[1],V0,P0, xaxis=:log, yaxis=:log,
        xlabel=L"\mathrm{Log}(V/V_0)", ylabel=L"\mathrm{Log}(P/P_0)",
        xlim=(0.04,1.1), ylim=(0.8, maximum(P0)*1.1),
        size=(720,540), legend=:none)
    plot!(plt[2],V0,P0,
        xlabel=L"V/V_0", ylabel=L"P/P_0",
        xlim=(0.0,1.1), ylim=(0, maximum(P0)*1.1),
        size=(720,540), legend=:none)
    plot(plt..., size=(1080,540))
end


function Vfun1(θ, CA_open, CA_close)
    #if abs(θ) .> π
    #    return 1
    #else
    Vcyl=(1 .+ 0.5*(rc-1) .* (R .+ 1.0 .- cos.(θ) .- ((R^2 .- (sin.(θ)).^2).^(1/2))))*Vc
    if (CA_open<θ) && (CA_close>θ)
        return (Vcyl+Vd+Vc)/(Vc+Vd);
    else
        return Vcyl/(Vc+Vd);
    end
    #end
end


###########################
# External Reactor

θ_open=(-40*π/180); θ_close=(40*π/180)
V1=Vfun1.(CA_rad, θ_open, θ_close)
P1=V1.^(-γ)
inds=1:length(P1)
rng_open=((CA_rad).>θ_open) .&& ((CA_rad).<θ_close)
inds_open=inds[rng_open]
P_open=P1[inds_open[1]-1]
P1[rng_open] .= P_open .* (V1[rng_open]./V1[inds_open[1]]).^(-γ)
P_close=P1[inds_open[end]]
P1_cut=copy(P1)
P1_cut[inds_open[1]]=NaN
P1_cut[inds_open[end]]=NaN


P1[abs.(CA_rad) .>= π] .= 1.0
P1_cut[abs.(CA_rad) .>= π] .= 1.0
using Plots, LaTeXStrings, Measures
begin
    plt=[Plots.plot() for _ in 1:2]

    plot!(plt[1],V1,P1_cut, xaxis=:log, yaxis=:log,
        xlabel=L"\mathrm{Log}(V/V_0)", ylabel=L"\mathrm{Log}(P/P_0)",
        ylim=(0.8, maximum(P1)*1.1),
        size=(720,540), legend=:right, label="Total V")
    plot!(plt[1],V0,P1,label="Cylinder V")

    plot!(plt[2],V1,P1_cut,
        xlabel=L"V/V_0", ylabel=L"P/P_0",
        ylim=(0.8, maximum(P1)*1.1),
        size=(720,540), legend=:right, label="Total V")
    plot!(plt[2],V0,P1,label="Cylinder V")
    plot(plt..., size=(1080,540), margin=10mm)
end


##########################
#Classic Otto

P2 = copy(P0)
inds_burned=inds[CA_rad .>= 0.0]
P2[inds_burned[1]]+=50;
P2[inds_burned]=P2[inds_burned[1]]*(V0[inds_burned]/V0[inds_burned[1]]).^(-γ)
P2[abs.(CA_rad) .>= π] .= 1.0
using Plots, LaTeXStrings
begin
    plt=[Plots.plot() for _ in 1:2]
    plot!(plt[1],V0,P2, xaxis=:log, yaxis=:log,
        xlabel=L"\mathrm{Log}(V/V_0)", ylabel=L"\mathrm{Log}(P/P_0)",
        xlim=(0.04,1.1), ylim=(0.8, maximum(P2)*1.1),
        size=(720,540), legend=:none)
    plot!(plt[2],V0,P2,
        xlabel=L"V/V_0", ylabel=L"P/P_0",
        xlim=(0.0,1.1), ylim=(0, maximum(P2)*1.1),
        size=(720,540), legend=:none)
    plot(plt..., size=(1080,540))
end

##########################
#Reactor with pressure rise/fall

for i in 1:2
    if i==1
        θ_open=(-40*π/180); θ_close=(30*π/180)
    else
        θ_open=(-30*π/180); θ_close=(40*π/180)
    end
    V3=Vfun1.(CA_rad, θ_open, θ_close)
    P3=V3.^(-γ)
    inds=1:length(P3)
    rng_open=((CA_rad).>θ_open) .&& ((CA_rad).<θ_close)
    inds_open=inds[rng_open]
    P_open=P3[inds_open[1]-1]
    P3[rng_open] .= P_open .* (V3[rng_open]./V3[inds_open[1]]).^(-γ)
    P_close=P3[inds_open[end]]
    inds_exp=inds[inds_open[end]+1:end]
    P3[inds_exp]=P3[inds_exp[1]-1]*(V0[inds_exp]/V0[inds_exp[1]-1]).^(-γ)



    P3[abs.(CA_rad) .>= π] .= 1.0
    begin
        plt=[Plots.plot() for _ in 1:2]

        plot!(plt[1],V0,P3, xaxis=:log, yaxis=:log,
            xlabel=L"\mathrm{Log}(V/V_0)", ylabel=L"\mathrm{Log}(P/P_0)",
            ylim=(0.5, maximum(P3)*1.1),
            size=(720,540), legend=:none)

        plot!(plt[2],V0,P3,
            xlabel=L"V/V_0", ylabel=L"P/P_0",
            ylim=(0, maximum(P3)*1.1),
            size=(720,540), legend=:none)
        plot(plt..., size=(1080,540), margin=10mm)
    end
end

###################################
# Pressure drop reactor with pumping work
θ_open=(-40*π/180); θ_close=(30*π/180)

V3=Vfun1.(CA_rad, θ_open, θ_close)
P3=V3.^(-γ)
inds=1:length(P3)
rng_open=((CA_rad).>θ_open) .&& ((CA_rad).<θ_close)
inds_open=inds[rng_open]
P_open=P3[inds_open[1]-1]
P3[rng_open] .= P_open .* (V3[rng_open]./V3[inds_open[1]]).^(-γ)
P_close=P3[inds_open[end]]
inds_exp=inds[inds_open[end]+1:end]
P3[inds_exp]=P3[inds_exp[1]-1]*(V0[inds_exp]/V0[inds_exp[1]-1]).^(-γ)



P3[CA_rad .<= -π] .= 1.0
P3[CA_rad .>= π] .= 2.0
V3=[V0;V0[1]]
P3=[P3;P3[1]]
begin
    plt=[Plots.plot() for _ in 1:2]

    plot!(plt[1],V3,P3, xaxis=:log, yaxis=:log,
        xlabel=L"\mathrm{Log}(V/V_0)", ylabel=L"\mathrm{Log}(P/P_0)",
        ylim=(0.5, maximum(P3)*1.1),
        size=(720,540), legend=:none)

    plot!(plt[2],V3,P3,
        xlabel=L"V/V_0", ylabel=L"P/P_0",
        ylim=(0, maximum(P3)*1.1),
        size=(720,540), legend=:none)
    plot(plt..., size=(1080,540), margin=10mm)
end
