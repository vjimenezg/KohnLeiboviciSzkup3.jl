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
            m=merge((w = Prices[2,1],
            Yc = Prices[1,1],
            ξ = Prices[3,1],
            Pk = Prices[4,1],
            Yk = Prices[5,1]),m)
        end
        solver=(x0 = [m.w m.Yc m.ξ m.Pk m.Yk],)

    else

        if s.load_prices ==0
            m=merge((w = 0.035690209553739,
            Φ_h = 0.16403,
            ξ =0.369060692554216,
            Pk = 0.975043532756500),m)
        elseif s.load_prices ==1
            Prices=readdlm("KLS3_prices.txt",'\t',Float64,'\n')
            m=merge((w = Prices[2,1],
            Φ_h = Prices[1,1],
            ξ = Prices[3,1],
            Pk = Prices[4,1]),m)
        end
        solver=(x0 = [m.w m.Φ_h m.ξ m.Pk],)

    end

    if s.GE == 1
        results_GE =
         nlsolve((F,x) ->KLS3_GE_par!(F,x,m,s,r),log.(solver.x0),autodiff = :forward, method=s.method, xtol=s.xtol_GE, ftol=s.ftol_GE, s.show_trace_GE)
        solver=merge((z=results_GE.zero,),solver)
        mcc_0, m_0, r_0, s_0, sim_0 = KLS3_GE_par(solver.z,m,s,r)
        solver=merge((mcc_0=mcc_0,),solver)
    elseif ge.GE ==0
        mcc_0, m_0, r_0, s_0, sim_0 =KLS3_GE_par(log.(solver.x0),m,s,r)
        solver=merge((mcc_0=mcc_0,),solver)
    end



    ## Step 2: Solve final steady state
     if s.transition == 1

        #Shocks
        m=merge((r = m.rv[end],
        Pf = m.Pfv[end],
        Yf = m.Yfv[end],
        β=m.β_v[end],
        θ = m.θ_v[end],
        τ = m.τ_v[end],
        δ = m.δ_v[end],
        τ_x = m.τ_x_v[end],
        τ_m_c = m.τ_m_c_v[end],
        τ_m_k = m.τ_m_k_v[end]),m)

         if s.tariffsincome == 1

            if s.load_prices ==0  &&  s.PE==0
                 m=merge((w = 0.021645899110592,
                Yc = 0.075394893630123,
                ξ = 0.280906535178547,
                Pk = 0.969892818482190,
                Yk = 0.065988921968899),m)
            elseif s.load_prices ==1
                Prices=readdlm("KLS3_prices.txt",'\t',Float64,'\n')
                m=merge((w = Prices[2,end],
                Yc = Prices[1,end],
                ξ = Prices[3,end],
                Pk = Prices[4,end],
                Yk = Prices[5,end]),m)
            elseif s.load_prices ==0  &&  s.PE==1
                m=merge((w = m_0.w,
                Yc = m_0.Yc,
                ξ = m_0.ξ,
                Pk = m_0.Pk,
                Yk = m_0.Yk),m)
            end
            solver=(x0 = [m.w m.Yc m.ξ m.Pk m.Yk],)

        else

            if s.load_prices ==0  &&  s.PE==0
                m=merge((w = 0.035690209553739,
                ϕ_h = 0.16403,
                ξ =0.369060692554216,
                Pk = 0.975043532756500),m)
            elseif s.load_prices ==1
                Prices=readdlm("KLS3_prices.txt",'\t',Float64,'\n')
                m=merge((w = Prices[2,1],
                ϕ_h = Prices[1,1],
                ξ = Prices[3,1],
                Pk = Prices[4,1]),m)
            elseif s.load_prices ==0 && s.PE==0
                m=merge((w = m_0.w,
                ϕ_h = m_0.ϕ_h,
                ξ = m_0.ξ,
                Pk = m_0.Pk),m)
            end

            solver=(x0 = [m.w m.Φ_h m.ξ m.Pk],)
        end
        if s.GE_end == 1  && s.PE == 0
            results_GE =
             nlsolve((F,x) ->KLS3_GE_par!(F,x,m,s,r),log.(solver.x0),autodiff = :forward, method=s.method, xtol=s.xtol_GE, ftol=s.ftol_GE, s.show_trace_GE)
            solver=merge((z=results_GE.zero,),solver)
            mcc_end, m_end, r_end, s_end, sim_end = KLS3_GE_par(solver.z,m,s,r)
            solver=merge((mcc_end=mcc_end,),solver)
        else
            mcc_end, m_end, r_end, s_end, sim_end =KLS3_GE_par(log.(solver.x0),m,s,r)
            solver=merge((mcc_end=mcc_end,),solver)
        end

        # Real variables
        sim_end=merge((X_Laspeyres = sum(sim_end.measure.*m_0.ξ.*r_0.pf.*r_end.yf),
        D_Laspeyres = sum(sim_end.measure.* r_0.pd.*r_end.yd),
        Sales_Laspeyres = sum(sim_end.measure.*(r_0.pd.*r_end.yd + m_0.ξ.*r_0.pf.*r_end.yf)),
        GDP_Laspeyres = sum(sim_end.measure.*(r_0.pd.*r_end.yd + m_0.ξ*r_0.pf.*r_end.yf -r_end.m*m_0.Pk))),sim_end)


end
