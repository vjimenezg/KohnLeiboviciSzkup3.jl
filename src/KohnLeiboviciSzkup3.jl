module KohnLeiboviciSzkup3

## ##################################################
# This code solves the model in
# Financial Development and Trade Liberalization
# Kohn, Leibovici, Szkup 2020
####################################################

clearconsole()

# Dependencies

using LinearAlgebra
using DifferentialEquations, Sundials, SimpleDifferentialOperators, DiffEqCallbacks
using Delimited Files, DataFrames, DataFramesMeta, CSV, CSVFiles, JSON # results caching
using Interpolations, QuadGK # integration
using NLsolve # root-finding
using NamedTupleTools, Parameters # named tuples
using Roots
using QuantEcon
using Base
# Model files

include("parameters_settings.jl")
include("utils.jl")

#Step 1: Solve initial steady state

    if t_inc.tariffsincome == 1
        if ad_opt.load_prices ==0
            p_ext=@with_kw (w= 0.021645899110592,
            Yc = 0.075394893630123,
            ξ = 0.280906535178547,
            Pk = 0.969892818482190,
            Yk = 0.065988921968899)
        elseif ad_opt.load_prices ==1
            Prices=readdlm("KLS3_prices.txt",'\t',Float64,'\n')
            p_ext=@with_kw (w = Prices(2,1)),
            Yc = Prices(1,1),
            ξ = Prices(3,1),
            Pk = Prices(4,1),
            Yk = Prices(5,1))
        end
        x0 = [p_ext.w p_ext.Yc p_ext.ξ p_ext.Pk p_ext.Yk];

    else

        if ad_opt.load_prices ==0
            p_ext=@with_kw (w = 0.035690209553739,
            Φ_h = 0.16403,
            ξ =0.369060692554216,
            Pk = 0.975043532756500)
        elseif ad_opt.load_prices ==1
            Prices=readdlm("KLS3_prices.txt",'\t',Float64,'\n')
            p_ext=@with_kw (w = Prices(2,1),
            Φ_h = Prices(1,1),
            ξ = Prices(3,1),
            Pk = Prices(4,1))
        end
        x0 = [p_ext.w p_ext.Φ_h p_ext.ξ p_ext.Pk];

    end
# Didn't recover fval nor exitflag, not sure if this is a problem [see https://github.com/JuliaNLSolvers/NLsolve.jl/blob/master/src/nlsolve/solver_state_results.jl]

    if ge.GE == 1
        results_GE =
         nlsolve((F,x) ->KLS3_GE_par(F,x,),log(x0),autodiff = :forward, method=op_GE.method, xtol=op_GE.xtol, ftol=op_GE.ftol, show_traceop_GE.show_trace)
        solver=@with_kw (x0=x0,
        z=results_GE.zero)

#From here onwards pending
        [solver.mcc_0, m_0, r_0, s_0, sim_0] = KLS3_GE_par(solver.z,m,s,r);
    elseif ge.GE ==0
        [solver.mcc_0, m_0, r_0, s_0, sim_0] = KLS3_GE_par(log(solver.x0),m,s,r);
    end

end
