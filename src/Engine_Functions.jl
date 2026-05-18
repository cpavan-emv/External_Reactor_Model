
function V_Vc(CA,CR,R)
    # this calculate V/V_clearance
    return 1 + 0.5 * (CR-1) * (R + 1 - cos(CA) - sqrt(R*R - sin(CA)*sin(CA)));
end
function dV_Vc(CA, CR,R)
    # this calculate d(V/V_clearance)/dCA
    return 0.5 * (CR-1) * (sin(CA) + 0.5/sqrt(R*R - sin(CA)*sin(CA))*2*sin(CA)*cos(CA));
end

function V_dVfunc(t, ω, CR,R)
    CA=ω*t + π
    V=V_Vc(CA,CR,R)
    dV=dV_Vc(CA,CR,R)*ω
    return V, dV
end

sigmoid(x, shift, scale)=(1 .+ exp.(-(x-shift)/scale)).^(-1)

function θ_dθ(t,ω, x0,s)
    # burn fraction function
    CA_deg=(ω*t-π)*180/π
    θ=sigmoid(CA_deg,x0,s)
    dθ=180*ω/π*(1/s)*θ.^2 .*exp.(-(CA_deg-x0)/s)
    return θ, dθ
end

