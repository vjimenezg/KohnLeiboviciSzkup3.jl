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
using NaNMath
# Model files

include("parameters_settings.jl")
include("utils.jl")
include("GE_functions.jl")

#Step 1: Solve initial steady state

    if s.tariffsincome == 1
        if s.load_prices ==0
            m=merge((w= 0.021645899110592,
            Yc = 0.075394893630123,
            ξ = 0.280906535178547,
            Pk = 0.969892818482190,
            Yk = 0.065988921968899),m)
        elseif s.load_prices ==1
            Prices=readdlm("KLS3_prices.txt",'\t',Float64,'\n')
            m=merge((w = Prices(2,1),
            Yc = Prices(1,1),
            ξ = Prices(3,1),
            Pk = Prices(4,1),
            Yk = Prices(5,1)),m)
        end
        x0 = [m.w m.Yc m.ξ m.Pk m.Yk]

    else

        if s.load_prices ==0
            m=merge((w = 0.035690209553739,
            Φ_h = 0.16403,
            ξ =0.369060692554216,
            Pk = 0.975043532756500),m)
        elseif s.load_prices ==1
            Prices=readdlm("KLS3_prices.txt",'\t',Float64,'\n')
            m=merge((w = Prices(2,1),
            Φ_h = Prices(1,1),
            ξ = Prices(3,1),
            Pk = Prices(4,1)),m)
        end
        solver=(x0 = [m.w m.Φ_h m.ξ m.Pk],)

    end

    if s.GE == 1
        results_GE =
         nlsolve((F,x) ->KLS3_GE_par(F,x,m,s,r),log(x0),autodiff = :forward, method=s.method, xtol=s.xtol_GE, ftol=s.ftol_GE, s.show_trace_GE)
        solver=merge((z=results_GE.zero,),solver)

#From here onwards pending
        [solver.mcc_0, m_0, r_0, s_0, sim_0] = KLS3_GE_par(solver.z,m,s,r);
    elseif ge.GE ==0
        [solver.mcc_0, m_0, r_0, s_0, sim_0] = KLS3_GE_par(log(solver.x0),m,s,r);
    end

end
