##############################
#=
The functions in this file are for describing the case of compression with 2 volumes
Volume A is the volume of fixed size and contains the catalyst (reactions allowed)
Volume B is the volume of variable size and contains no catalyst (no reactions)
The two volumes evolve separately but are solved in the same vector
This allows for a callback function to couple them
    (open a check valve when pressures are equal)
=#
##############################

function f_DC!(du, u, p, t)
    f_CV!((@view du[1:1+Nspec]), u[1:1+Nspec],p,t)
    f_Comp!((@view du[2+Nspec:end]), u[2+Nspec:end], p, t)
end

function evolve_DC(u0,param,tsave, condition=nothing)
    prob=ODE.ODEProblem(DE_Model.f_DC!,u0,(tsave[1],tsave[end]))
    if !isnothing(condition)
        cb=ODE.ContinuousCallback(condition, (integrator)->ODE.terminate!(integrator))
        sol = ODE.solve(prob,ODE.RadauIIA5(autodiff=false),
            p=param, reltol=1e-8, abstol=1e-10, saveat=tsave, callback=cb)
    else
        sol = ODE.solve(prob,ODE.RadauIIA5(autodiff=false),
            p=param, reltol=1e-8, abstol=1e-10, saveat=tsave)
    end
    y=permutedims(stack(sol.u))
    return sol.t, y
end

function Vars_Eq_condition(u,t,integrator, inds)
    return u[inds[1]] - u[inds[2]]
end


