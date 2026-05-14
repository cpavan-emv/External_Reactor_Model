##############################
#=
The functions in this file are for describing the case of compression with of a non-reacting gas
Evolving species conservation is unnecessary, but may make tracking species easier
Also allows re-using mass functions from CV reactor
=#
##############################
# the mass matrix functions are identical to the CV reactor versions
# only difference is the RHS
# the @inline is to attempt to copy the text exactly at compile time (hopefully reducing overhead)

function M_Comp(gasB)
    return @inline M_CV(gasB)
end

# Matrix D-F are all the same as the "common" version
function D_Comp(gasB)
    return @inline D_CV(gasB)
end

function B_Comp(gasB)
    return @inline B_CV(gasB)
end

function F_Comp(gasB)
    return @inline F_CV(gasB)
end

function setMassMat_Comp(u,p,t)
    gasB=p[2] # for passing to fluid props
    fluid_props.setTPX(gasB, (u[2], u[1], u[3:end]), spec_ind)
    return M_Comp(gasB)*
         F_Comp(gasB)*
         B_Comp(gasB)*
         D_Comp(gasB)
end

function f_Comp!(du, u, p, t)
    # this is the entire RHS of the function
    # should only be used with stiff equation solvers
    Mass=setMassMat_Comp(u,p,t)
    gasB=p[2] # for passing to fluid props
    V, dV = p[3](t) # function of t returning V and dV

    du .= Mass\(-dV*fluid_props.rho(gasB)/V *vcat([1.0; fluid_props.enthalpy(gasB)],
        fluid_props.Y(gasB)[spec_ind[1:end-1]]))

    return nothing
end
