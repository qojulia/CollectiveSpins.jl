import OrdinaryDiffEqCore, OrdinaryDiffEqLowOrderRK, DiffEqCallbacks

"""
integrate()
"""
function integrate(T::Vector, f:: Function, state0::S, fout::Function;
                    alg = OrdinaryDiffEqLowOrderRK.DP5(),
                    callback = nothing, kwargs...) where S

    if isa(state0, Vector{<:Real})
        x0 = state0
        N = length(state0)
        fout_diff = (u, t, integrator) -> fout(t, deepcopy(u))
    else
        x0 = state0.data
        N = state0.N
        fout_diff = (u, t, integrator) -> fout(t, S(N, deepcopy(u)))
    end

    out_type = pure_inference(fout, Tuple{eltype(T),typeof(state0)})
    out = DiffEqCallbacks.SavedValues(eltype(T),out_type)
    scb = DiffEqCallbacks.SavingCallback(fout_diff,out,saveat=T,
                                        save_everystep=false,
                                        save_start = false)

    prob = OrdinaryDiffEqCore.ODEProblem(f, x0, (T[1], T[end]))

    full_cb = OrdinaryDiffEqCore.CallbackSet(callback, scb)

    sol = OrdinaryDiffEqCore.solve(prob, alg;
            reltol=1.0e-6,
            abstol=1.0e-8,
            save_everystep = false,
            save_start = false,
            save_end = false,
            callback=full_cb,
            kwargs...)

    out.t, out.saveval

end

Base.@pure pure_inference(fout,T) = Core.Compiler.return_type(fout, T)
