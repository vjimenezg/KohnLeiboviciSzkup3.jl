module KohnLeiboviciSzkup3

## ##################################################
# This code solves the model in
# No Credit, No Gain: Trade Liberalization Dynamics, Production Inputs, and Financial Development
# Kohn, Leibovici, Szkup 2020
####################################################
clearconsole()

# Dependencies (commented the ones I'm not sure I need yet)

using LinearAlgebra
using Plots
using Distributions
using NLsolve # root-finding
using NamedTupleTools # named tuples
using QuantEcon #used to generate a markov chain
using Base
using Dates #used to save the date of different workspaces
using JLD2 #allows to save the workspace or individual variables in .jld2 format (akin to matlab's .mat)
using MAT #allows to save the workspace or individual variables in .mat format (useful for matlab comparisons and the welfare graphs)
using NaNMath
using Profile #used for benchmarking purposes
using OwnTime #used for benchmarking purposes
using Threads
# using Delimited Files, DataFrames, DataFramesMeta, CSV, CSVFiles, JSON # results caching
#using TimerOutputs.jl # allows to print nicely formatted tables with the timing of different sections of the code
#using GMT #3d histograms (not available yet)

# Model files
include("tauchen.jl")
include("parameters_settings.jl")
include("functions.jl")
include("functions_trans.jl")

say(what) = run(`osascript -e "say \"$(what)\""`, wait=false)
#function that lets you know when code finished running

#Step 1: Solve initial steady state

    if s.tariffsincome == 1
        if s.load_prices ==0
            m=merge(m,(w= 0.021645899110592,
            Yc = 0.075394893630123,
            ξ = 0.280906535178547,
            Pk = 0.969892818482190,
            Yk = 0.065988921968899))
        elseif s.load_prices ==1
            @load "KLS3_Prices.jld2" Prices
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
            ϕ_h = 0.16403,
            ξ =0.369060692554216,
            Pk = 0.975043532756500))
        elseif s.load_prices ==1
            @load "KLS3_Prices.jld2" Prices
            m=merge(m,(w = Prices[2,1],
            ϕ_h = Prices[1,1],
            ξ = Prices[3,1],
            Pk = Prices[4,1]))
        end
        solver=(x0 = [m.w m.ϕ_h m.ξ m.Pk],)

    end

@time begin
    if s.GE == 1
        results_GE =
         nlsolve((F,x) ->KLS3_GE_par!(F,x,m,s,r),log.(solver.x0),autodiff = :forward, method=s.method_GE, xtol=s.xtol_GE, ftol=s.ftol_GE, show_trace=s.show_trace_GE)
        solver=merge(solver,(z=results_GE.zero,))
        mcc_0, m_0, r_0, s_0, sim_0 = KLS3_GE_par(solver.z,m,s,r)
        solver=merge(solver,(mcc_0=mcc_0,))
    elseif ge.GE ==0
        mcc_0, m_0, r_0, s_0, sim_0 =KLS3_GE_par(log.(solver.x0),m,s,r)
        solver=merge(solver,(mcc_0=mcc_0,))

         #Testing with MATLAB solver prices
        #mcc_0_matlabprices, m_0_matlabprices, r_0_matlabprices, s_0_matlabprices, sim_0_matlabprices =KLS3_GE_par(log.([Prices_MATLAB[2,1] Prices_MATLAB[1,1] Prices_MATLAB[3,1] Prices_MATLAB[4,1] Prices_MATLAB[5,1]]),m,s,r)
    end
say("Initial Steady State Finished - you'd better come and take a look....")
end


    if s.tariffsincome_PE_flag==1
        s=merge(s,(tariffsincome_PE = 1,))
        m=merge(m,(tariffsincome = m_0.tariffsincome,))
    end

# Step 2: Solve final steady state
    if s.transition == 1

        #Shocks
        m=merge(m,(r = m.rv[end],
        Pf = m.Pfv[end],
        Yf = m.Yfv[end],
        β = m.β_v[end],
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
                @load "KLS3_Prices.jld2" Prices
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
                @load "KLS3_Prices.jld2" Prices
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

            solver=(x0 = [m.w m.ϕ_h m.ξ m.Pk],)
        end
        @time begin
            if s.GE_end == 1  && s.PE == 0
                results_GE = nlsolve((F,x) ->KLS3_GE_par!(F,x,m,s,r),log.(solver.x0), autodiff = :forward, method=s.method_GE, xtol=s.xtol_GE, ftol=s.ftol_GE, show_trace=s.show_trace_GE)
                solver=merge(solver,(z=results_GE.zero,))
                mcc_end, m_end, r_end, s_end, sim_end = KLS3_GE_par(solver.z,m,s,r)
                solver=merge(solver,(mcc_end=mcc_end,))
            else
            mcc_end, m_end, r_end, s_end, sim_end =KLS3_GE_par(log.(solver.x0),m,s,r)
            solver=merge(solver,(mcc_end=mcc_end,))
            end
            say("Final Steady State Finished - you'd better come and take a look....")
        end

            #Saving prices
            #file=matopen("Prices_Julia_2ndSS.mat","w")
            #write(file,"Prices_Julia_2SS",solver.z)
            #close(file)

            #Testing with MATLAB solver prices
            #file = matopen("Prices_MATLAB_2ndSS.mat")
            #Prices_MATLAB_2SS=read(file, "Prices_MATLAB_2SS")
            #close(file)
            #mcc_end_matlab, m_end_matlab, r_0_matlab, s_end_matlab, sim_end_matlab =KLS3_GE_par(Prices_MATLAB_2SS,m,s,r)


        # Real variables
        sim_end=merge(sim_end,(X_Laspeyres = sum(sim_end.measure.*m_0.ξ.*r_0.pf.*r_end.yf),
        D_Laspeyres = sum(sim_end.measure.* r_0.pd.*r_end.yd),
        Sales_Laspeyres = sum(sim_end.measure.*(r_0.pd.*r_end.yd + m_0.ξ.*r_0.pf.*r_end.yf)),
        GDP_Laspeyres = sum(sim_end.measure.*(r_0.pd.*r_end.yd + m_0.ξ*r_0.pf.*r_end.yf -r_end.m*m_0.Pk))))

 # STEP 3: guess sequence of aggregate prices for N period, with N sufficiently large (this was in a separate script in matlab, KLS3_trans_obj)

 # Store results from initial and final steady states

rt=Vector{Any}(undef,s.N)

rt[1] = merge(r_0,(measure=sim_0.measure,))
rt[s.N] = merge(r_end,(measure=sim_end.measure,))
 #r_0=nothing
 #r_end=nothing

 # Guess sequence of aggregate prices for N period, with N sufficiently large
 #ϕht = zeros(s.N,1)
 m=merge(m,(Pkt = zeros(s.N,1),
 wt = zeros(s.N,1),
 ξt = zeros(s.N,1),
 Ykt = zeros(s.N,1),
 Yct = zeros(s.N,1)))

    if s.transition_AD == 1
    m=merge(m,(tariffsincomet = zeros(Real,s.N,1),))
    else
    m=merge(m,(tariffsincomet = zeros(s.N,1),))
    end

    if s.tariffsincome == 1
     # Guess
     #m.ϕht[1] = m_0.ϕ_h
     m.wt[1] = m_0.w
     m.ξt[1] = m_0.ξ
     m.Pkt[1] = m_0.Pk
     m.Yct[1] = m_0.Yc
     m.Ykt[1] = m_0.Yk

        if s.PE == 0
         #m.ϕht[s.N]=m_end.ϕ_h
         m.wt[s.N]=m_end.w
         m.ξt[s.N]=m_end.ξ
         m.Pkt[s.N]=m_end.Pk
         m.Yct[s.N] = m_end.Yc
         m.Ykt[s.N] = m_end.Yk

        elseif s.PE ==1

         #m.ϕht[s.N]=m_0.ϕ_h
         m.wt[s.N]=m_0.w
         m.ξt[s.N]=m_0.ξ
         m.Pkt[s.N]=m_end.Pk #This was in the MATLAB code, is it a typo and should it be m_0.Pk as the rest of the prices?
         m.Yct[s.N] = m_0.Yc
         m.Ykt[s.N] = m_0.Yk

        end

     period2_change = 0.5

     init =m.wt[1]+ (m.wt[s.N]-m.wt[1])*period2_change
     expo =log.(m.wt[s.N]./init)./log.(s.N)   #m.wt[1]*(s.N)^expo=m.wt[s.N] =>expo=log.(m.wt[s.N]/m.wt[1])/log.(s.N)
     m.wt[2:s.N-1].=init*(2:s.N-1).^expo

     init =m.Yct[1]+ (m.Yct[s.N]-m.Yct[1])*period2_change
     expo =log.(m.Yct[s.N]./init)./log.(s.N) #m.wt[1]*(s.N)^expo=m.wt[s.N] =>expo=log.(m.wt[s.N]/m.wt[1])/log.(s.N)
     m.Yct[2:s.N-1].=init*(2:s.N-1).^expo

     init =m.Ykt[1]+ (m.Ykt[s.N]-m.Ykt[1])*period2_change
     expo =log.(m.Ykt[s.N]./init)./log.(s.N)   #m.wt[1]*(s.N).^expo=m.wt[s.N] =>expo=log.(m.wt[s.N]./m.wt[1])/log.(s.N)
     m.Ykt[2:s.N-1].=init*(2:s.N-1).^expo

     m.Pkt[2:s.N-1].=m_end.Pk
     m.ξt[2:s.N-1].=m_end.ξ

    Guess = [m.Yct[2:s.N-1] m.wt[2:s.N-1] m.ξt[2:s.N-1] m.Pkt[2:s.N-1] m.Ykt[2:s.N-1]]

        if s.load_prices == 1
        @load "KLS3_Prices.jld2" Prices
         Guess = [Prices[1,2:s.N-1] Prices[2,2:s.N-1] Prices[3,2:s.N-1] Prices[4,2:s.N-1] Prices[5,2:s.N-1]]
        end
    end

 ## STEP 4: Solving for the GE prices and quantities
@time begin
        if s.transition_GE == 1
            Guess=log.(Guess)
            if s.transition_AD == 1
                results_trans = nlsolve((F,x) -> KLS3_transition_vec2!(F,x,m,r,s,rt),Guess, autodiff = :forward, method=s.method_trans, factor=0.0001, autoscale=false, xtol=s.xtol_trans, ftol=s.ftol_trans,iterations=s.MaxIter_trans, show_trace=s.show_trace_trans,store_trace=true)
            else
                results_trans = nlsolve((F,x) -> KLS3_transition_vec2!(F,x,m,r,s,rt),Guess, method=s.method_trans, factor=1, autoscale=true, xtol=s.xtol_trans, ftol=s.ftol_trans,iterations=40, show_trace=s.show_trace_trans,store_trace=true)
            end
            Prices_sol=results_trans.zero
            mc, m, r, sim_fun, rt = KLS3_transition_vec2(Prices_sol,m,r,s,rt)
        else
            Prices_sol = log.(Guess)
            mc, m, r, sim_fun, rt = KLS3_transition_vec2(Prices_sol,m,r,s,rt)
        end
    say("Transition Finished - you'd better come and take a look....")
    end
end

#Trying with MATLAB prices
#file = matopen("Prices_Transition_MATLAB_Julia2ndSS.mat")
#Prices_sol_MATLAB=read(file, "Prices_sol")
#close(file)
#results_trans_matlab = nlsolve((F,x) -> KLS3_transition_vec2!(F,x,m,r,s,rt),Prices_sol_MATLAB, method=s.method_trans, factor=0.001, autoscale=false, xtol=s.xtol_trans, ftol=s.ftol_trans,iterations=40, show_trace=s.show_trace_trans,store_trace=true)
#Prices_sol_MATLABGuess=results_trans_matlab.zero
#mc_matlab, m_matlab, r_matlab, sim_fun_matlab, rt_matlab = KLS3_transition_vec2(Prices_sol_MATLABGuess,m,r,s,rt)


 ## Welfare analysis
 if s.welfare==1 #FIGURES PENDING
include("welfare.jl")
 end

 ## End

 if s.save_prices==1
     if s.tariffsincome == 0
         Prices[1,:] = [m.ϕht[1] exp.(Prices_sol[1:s.N-2])' m.ϕht[s.N]]
         Prices[2,:] = [m.wt[1] exp.(Prices_sol[s.N-1:2*s.N-4])' m.wt[s.N]]
         Prices[3,:] = [m.ξt[1] exp.(Prices_sol[2*s.N-3:3*s.N-6])' m.ξt[s.N]]
         Prices[4,:] = [m.Pkt[1] exp.(Prices_sol[3*s.N-5:4*s.N-8])' m.Pkt[s.N]]
     elseif s.tariffsincome == 1
         Prices[1,:] = [m.Yct[1] exp.(Prices_sol[1:s.N-2])' m.Yct[s.N]]
         Prices[2,:] = [m.wt[1] exp.(Prices_sol[s.N-1:2*s.N-4])' m.wt[s.N]]
         Prices[3,:] = [m.ξt[1] exp.(Prices_sol[2*s.N-3:3*s.N-6])' m.ξt[s.N]]
         Prices[4,:] = [m.Pkt[1] exp.(Prices_sol[3*s.N-5:4*s.N-8])' m.Pkt[s.N]]
         Prices[5,:] = [m.Ykt[1] exp.(Prices_sol[4*s.N-7:end])' m.Ykt[s.N]]
     end
    @save "KLS3_prices.jld2" Prices
 end

 if s.save_workspace==1
     rt = nothing
     r_0 = nothing
     r_end = nothing
     @save "workspace_$(Dates.year(now()))_$(Dates.monthname(now()))_$(Dates.day(now()))_$(Dates.hour(now()))_$(Dates.minute(now()))"

 end


end
