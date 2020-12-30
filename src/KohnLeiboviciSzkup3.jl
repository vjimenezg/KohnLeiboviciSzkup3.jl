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
include("functions.jl")

#Step 1: Solve initial steady state

    if s.tariffsincome == 1
        if s.load_prices ==0
            m=merge(m,(w= 0.021645899110592,
            Yc = 0.075394893630123,
            ξ = 0.280906535178547,
            Pk = 0.969892818482190,
            Yk = 0.065988921968899))
        elseif s.load_prices ==1
            Prices=readdlm("KLS3_prices.txt",'\t',Float64,'\n')
            m=merge(m,(w = Prices[2,1],
            Yc = Prices[1,1],
            ξ = Prices[3,1],
            Pk = Prices[4,1],
            Yk = Prices[5,1]))
        end
        solver=(x0 = [m.w m.Yc m.ξ m.Pk m.Yk],)

    else

        if s.load_prices ==0
            m=merge(m,(w = 0.035690209553739,
            Φ_h = 0.16403,
            ξ =0.369060692554216,
            Pk = 0.975043532756500))
        elseif s.load_prices ==1
            Prices=readdlm("KLS3_prices.txt",'\t',Float64,'\n')
            m=merge(m,(w = Prices[2,1],
            Φ_h = Prices[1,1],
            ξ = Prices[3,1],
            Pk = Prices[4,1]))
        end
        solver=(x0 = [m.w m.Φ_h m.ξ m.Pk],)

    end

    if s.GE == 1
        results_GE =
         nlsolve((F,x) ->KLS3_GE_par!(F,x,m,s,r),log.(solver.x0),autodiff = :forward, method=s.method, xtol=s.xtol_GE, ftol=s.ftol_GE, s.show_trace_GE)
        solver=merge(solver,(z=results_GE.zero,))
        mcc_0, m_0, r_0, s_0, sim_0 = KLS3_GE_par(solver.z,m,s,r)
        solver=merge(solver,(mcc_0=mcc_0,))
    elseif ge.GE ==0
        mcc_0, m_0, r_0, s_0, sim_0 =KLS3_GE_par(log.(solver.x0),m,s,r)
        solver=merge(solver,(mcc_0=mcc_0,))
    end



    ## Step 2: Solve final steady state
     if s.transition == 1

        #Shocks
        m=merge(m,(r = m.rv[end],
        Pf = m.Pfv[end],
        Yf = m.Yfv[end],
        β=m.β_v[end],
        θ = m.θ_v[end],
        τ = m.τ_v[end],
        δ = m.δ_v[end],
        τ_x = m.τ_x_v[end],
        τ_m_c = m.τ_m_c_v[end],
        τ_m_k = m.τ_m_k_v[end]))

         if s.tariffsincome == 1

            if s.load_prices ==0  &&  s.PE==0
                 m=merge(m,(w = 0.021645899110592,
                Yc = 0.075394893630123,
                ξ = 0.280906535178547,
                Pk = 0.969892818482190,
                Yk = 0.065988921968899))
            elseif s.load_prices ==1
                Prices=readdlm("KLS3_prices.txt",'\t',Float64,'\n')
                m=merge(m,(w = Prices[2,end],
                Yc = Prices[1,end],
                ξ = Prices[3,end],
                Pk = Prices[4,end],
                Yk = Prices[5,end]))
            elseif s.load_prices ==0  &&  s.PE==1
                m=merge(m,(w = m_0.w,
                Yc = m_0.Yc,
                ξ = m_0.ξ,
                Pk = m_0.Pk,
                Yk = m_0.Yk))
            end
            solver=(x0 = [m.w m.Yc m.ξ m.Pk m.Yk],)

        else

            if s.load_prices ==0  &&  s.PE==0
                m=merge(m,(w = 0.035690209553739,
                ϕ_h = 0.16403,
                ξ =0.369060692554216,
                Pk = 0.975043532756500))
            elseif s.load_prices ==1
                Prices=readdlm("KLS3_prices.txt",'\t',Float64,'\n')
                m=merge(m,(w = Prices[2,1],
                ϕ_h = Prices[1,1],
                ξ = Prices[3,1],
                Pk = Prices[4,1]))
            elseif s.load_prices ==0 && s.PE==0
                m=merge(m,(w = m_0.w,
                ϕ_h = m_0.ϕ_h,
                ξ = m_0.ξ,
                Pk = m_0.Pk))
            end

            solver=(x0 = [m.w m.Φ_h m.ξ m.Pk],)
        end
        if s.GE_end == 1  && s.PE == 0
            results_GE =
             nlsolve((F,x) ->KLS3_GE_par!(F,x,m,s,r),log.(solver.x0),autodiff = :forward, method=s.method, xtol=s.xtol_GE, ftol=s.ftol_GE, s.show_trace_GE)
            solver=merge(solver,(z=results_GE.zero,))
            mcc_end, m_end, r_end, s_end, sim_end = KLS3_GE_par(solver.z,m,s,r)
            solver=merge(solver,(mcc_end=mcc_end,))
        else
            mcc_end, m_end, r_end, s_end, sim_end =KLS3_GE_par(log.(solver.x0),m,s,r)
            solver=merge(solver,(mcc_end=mcc_end,))
        end

        # Real variables
        sim_end=merge(sim_end,(X_Laspeyres = sum(sim_end.measure.*m_0.ξ.*r_0.pf.*r_end.yf),
        D_Laspeyres = sum(sim_end.measure.* r_0.pd.*r_end.yd),
        Sales_Laspeyres = sum(sim_end.measure.*(r_0.pd.*r_end.yd + m_0.ξ.*r_0.pf.*r_end.yf)),
        GDP_Laspeyres = sum(sim_end.measure.*(r_0.pd.*r_end.yd + m_0.ξ*r_0.pf.*r_end.yf -r_end.m*m_0.Pk))))

 ## STEP 3: guess sequence of aggregate prices for N period, with N sufficiently large (this was in a separate script in matlab, KL3_trans_obj)

 # Store results from initial and final steady states [PENDING]

 rt{1,1} = r_0
 rt{1,1}.measure= sim_0.measure
 rt{s.N,1} = r_end
 rt{s.N,1}.measure= sim_end.measure
 r_0=nothing
 r_end=nothing

 # Guess sequence of aggregate prices for N period, with N sufficiently large
 #ϕht = zeros(s.N,1)
 Pkt = zeros(s.N,1)
 wt = zeros(s.N,1)
 ξt = zeros(s.N,1)
 Ykt = zeros(s.N,1)
 Yct = zeros(s.N,1)
 tariffsincomet= zeros(s.N,1)

 if s.tariffsincome == 1
     # Guess
     #ϕht[1] = m_0.ϕ_h
     wt[1] = m_0.w
     ξt[1] = m_0.ξ
     Pkt[1] = m_0.Pk
     Yct[1] = m_0.Yc
     Ykt[1] = m_0.Yk

     if s.PE == 0
         #ϕht[s.N]=m_end.ϕ_h
         wt[s.N]=m_end.w
         ξt[s.N]=m_end.ξ
         Pkt[s.N]=m_end.Pk
         Yct[s.N] = m_end.Yc
         Ykt[s.N] = m_end.Yk

     elseif s.PE ==1

         #ϕht[s.N]=m_0.ϕ_h
         wt[s.N]=m_0.w
         ξt[s.N]=m_0.ξ
         Pkt[s.N]=m_end.Pk
         Yct[s.N] = m_0.Yc
         Ykt[s.N] = m_0.Yk

     end


     period2_change = 0.5

     init =wt(1)+ (wt[s.N]-wt[1])*period2_change
     expo =log.(wt[s.N]./init)./log.(s.N)   #wt[1]*(s.N)^expo=wt[s.N] =>expo=log.(wt[s.N]/wt[1])/log.(s.N)
     wt[2:s.N-1]=init*(2:s.N-1).^expo

     init =Yct[1]+ (Yct[s.N]-Yct[1])*period2_change
     expo =log.(Yct[s.N]./init)./log.(s.N) #wt[1]*(s.N)^expo=wt[s.N] =>expo=log.(wt[s.N]/wt[1])/log.(s.N)
     Yct[2:s.N-1]=init*(2:s.N-1).^expo

     init =Ykt[1]+ (Ykt[s.N]-Ykt[1])*period2_change
     expo =log.(Ykt[s.N]./init)./log.(s.N)   #wt[1]*(s.N).^expo=wt[s.N] =>expo=log.(wt[s.N]./wt[1])/log.(s.N)
     Ykt[2:s.N-1]=init*(2:s.N-1).^expo

     Pkt[2:s.N-1]=m_end.Pk
     ξt[2:s.N-1]=m_end.ξ

    Guess = [Yct[2:s.N-1] wt[2:s.N-1] ξt[2:s.N-1] Pkt[2:s.N-1] Ykt[2:s.N-1]]

    if s.load_prices == 1
        Prices=readdlm("KLS3_prices.txt",'\t',Float64,'\n')
         Guess = [Prices[1,2:s.N-1] Prices[2,2:s.N-1] Prices[3,2:s.N-1] Prices[4,2:s.N-1] Prices[5,2:s.N-1]]
     end

     m=merge(m,(Pkt=Pkt,wt=wt,ξt=ξt,Ykt=Ykt,Yct=Yct,tariffsincomet=tariffsincomet))
     #m=merge(m,(ϕht=ϕht,))
 end

 ## STEP 4: Solving for the GE prices and quantities

     if s.transition_GE== 1
         # guess
         Guess=log.(Guess)
         # solution
         solve_prices = @(Prices)KLS3_transition_vec2(Prices,m,r,s,rt); #Function that takes prices as inputs and market clearing values as output

         ### FROM HERE ONWARDS PENDING, WORKING ON KLS3_transition_vect2 ###
         [Prices_sol, mcc_sol, exit_sol] = fsolve(solve_prices,Guess,s.options_trans)
         [mc, m, r, sim_fun, rt] = KLS3_transition_vec2(Prices_sol,m,r,s,rt)
     else

         Guess=log.(Guess)
         Prices_sol = Guess
         mc, m, r, sim_fun, rt = KLS3_transition_vec2(Prices_sol,m,r,s,rt)

     end


 end


end
