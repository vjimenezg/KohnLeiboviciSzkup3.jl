module KohnLeiboviciSzkup3 #introduces a new global scope

## ##################################################
# This code solves the model in
# No Credit, No Gain: Trade Liberalization Dynamics, Production Inputs, and Financial Development
# Kohn, Leibovici, Szkup 2021
####################################################

# Dependencies (packages used)
using Base
using LinearAlgebra
using Distributions
using Plots # used to plot graphs
using Parameters # used to create named tuples
using NLsolve # (multivariable) root-finder [main solver]
using LeastSquaresOptim # alternative solver [transition only]
using QuantEcon #used to generate a markov chain
using Dates #used to save the date of different workspaces
using JLD2 #allows to save the workspace or individual variables in .jld2 format (akin to matlab's .mat)
using MAT #allows to save the workspace or individual variables in .mat format (useful for matlab comparisons and the welfare graphs)

### Benchmarking packages
#using Profile
#using TimerOutputs.jl # allows to print nicely formatted tables with the timing of different sections of the code

# Model files
include("parameters_settings.jl") #contains the default setup of settings and parameters (+tauchen function)
include("functions.jl") #contains all functions used in the computation of the SS
include("functions_trans.jl") #contains (additional) functions used in the transition

p=parameters_default() #instantiating the parameters with default values
s=settings_default() #instantiating the settings with default values

say(what) = run(`osascript -e "say \"$(what)\""`, wait=false)
#function that lets you know when code finished running

#Step 1: Solve initial steady state

    if s.tariffsincome == 1
        if s.load_prices ==0
            const g1=(w= 0.021645899110592, # In Julia, it's better perfomance-wise not to use global variables such as this. If they are used and they are constants (no not change through the code), labeling them as such using "const" helps mitigate the effect on performance.
            Yc = 0.075394893630123,
            ξ = 0.280906535178547,
            Pk = 0.969892818482190,
            Yk = 0.065988921968899)
        elseif s.load_prices ==1
            @load "KLS3_Prices.jld2" Prices
            const g1=(w = Prices[2,1],
            Yc = Prices[1,1],
            ξ = Prices[3,1],
            Pk = Prices[4,1],
            Yk = Prices[5,1])
        end
        x0 = [g1.w g1.Yc g1.ξ g1.Pk g1.Yk]

    else

        if s.load_prices ==0
            const g1=(w = 0.035690209553739,
            ϕ_h = 0.16403,
            ξ =0.369060692554216,
            Pk = 0.975043532756500)
        elseif s.load_prices ==1
            @load "KLS3_Prices.jld2" Prices
            const g1=(w = Prices[2,1],
            ϕ_h = Prices[1,1],
            ξ = Prices[3,1],
            Pk = Prices[4,1])
        end
        x0 = [g1.w g1.ϕ_h g1.ξ g1.Pk]

    end

@time begin
    if s.GE == 1
        results_GE =
         nlsolve((F,x) ->KLS3_GE_par!(F,x,p,s),log.(x0),autodiff = :forward, method=s.method_GE, xtol=s.xtol_GE, ftol=s.ftol_GE, show_trace=s.show_trace_GE)
         mcc_0, p_o_0, sim_0 = KLS3_GE_par(results_GE.zero,p,s)
    elseif s.GE ==0
        mcc_0, p_o_0, sim_0 =KLS3_GE_par(log.(x0),p,s)
    end
    say("Initial Steady State Finished - you'd better come and take a look....")
end

#clearing previously defined global vars (so they don't retain their values when computing the 2nd SS)
guessM=nothing
guessV=nothing

    if s.tariffsincome_PE_flag==1
        s=settings_default(tariffsincome_PE = 1) #note that here what you do is instantiate "setting_default" as before but changing the default value of tariffsincome_PE (instead of changing the value in the previous instantiation, which cannot be done)
    end

# Step 2: Solve final steady state
    if s.transition == 1

        #Shocks
        if s.tariffsincome_PE_flag==1
            p=parameters_default(tariffsincome=p_o_0.tariffsincome,
            #The change in p.tarrifsincome was originally in the previous block of code. However, the fact that what we do here is instantiate p again, it's convenient to do all the changes to p together (otherwise, we'd have to add the change in tariffsincome again)
            r = p.rv[end],
            Pf = p.Pfv[end],
            Yf = p.Yfv[end],
            β = p.β_v[end],
            θ = p.θ_v[end],
            τ = p.τ_v[end],
            δ = p.δ_v[end],
            τ_x = p.τ_x_v[end],
            τ_m_c = p.τ_m_c_v[end],
            τ_m_k = p.τ_m_k_v[end])
        else
            p=parameters_default(r = p.rv[end],
            Pf = p.Pfv[end],
            Yf = p.Yfv[end],
            β = p.β_v[end],
            θ = p.θ_v[end],
            τ = p.τ_v[end],
            δ = p.δ_v[end],
            τ_x = p.τ_x_v[end],
            τ_m_c = p.τ_m_c_v[end],
            τ_m_k = p.τ_m_k_v[end])
        end

         if s.tariffsincome == 1

            if s.load_prices ==0  &&  s.PE==0
                 const g2=(w = 0.021645899110592,
                Yc = 0.075394893630123,
                ξ = 0.280906535178547,
                Pk = 0.969892818482190,
                Yk = 0.065988921968899)
            elseif s.load_prices ==1
                @load "KLS3_Prices.jld2" Prices
                const g2=(w = Prices[2,end],
                Yc = Prices[1,end],
                ξ = Prices[3,end],
                Pk = Prices[4,end],
                Yk = Prices[5,end])
            elseif s.load_prices ==0  &&  s.PE==1
                const g2=(w = p_o_0.w,
                Yc = p_o_0.Yc,
                ξ = p_o_0.ξ,
                Pk = p_o_0.Pk,
                Yk = p_o_0.Yk)
            end
            x0 = [g2.w g2.Yc g2.ξ g2.Pk g2.Yk]

        else

            if s.load_prices ==0  &&  s.PE==0
                const g2=(w = 0.035690209553739,
                ϕ_h = 0.16403,
                ξ =0.369060692554216,
                Pk = 0.975043532756500)
            elseif s.load_prices ==1
                @load "KLS3_Prices.jld2" Prices
                const g2=(w = Prices[2,1],
                ϕ_h = Prices[1,1],
                ξ = Prices[3,1],
                Pk = Prices[4,1])
            elseif s.load_prices ==0 && s.PE==0
                const g2=(w = p_o_0.w,
                ϕ_h = p_o_0.ϕ_h,
                ξ = p_o_0.ξ,
                Pk = p_o_0.Pk)
            end

            x0 = [g2.w g2.ϕ_h g2.ξ g2.Pk]
        end

        @time begin
            if s.GE_end == 1  && s.PE == 0
                results_GE = nlsolve((F,x) ->KLS3_GE_par!(F,x,p,s),log.(x0), autodiff = :forward, method=s.method_GE, xtol=s.xtol_GE, ftol=s.ftol_GE, show_trace=s.show_trace_GE)
                mcc_end, p_o_end, sim_end = KLS3_GE_par(results_GE.zero,p,s)
            else
                mcc_end, p_o_end, sim_end =KLS3_GE_par(log.(x0),p,s)
            end
            say("Final Steady State Finished - you'd better come and take a look....")
        end

        #clearing previously defined global vars
        guessV=nothing
        guessM=nothing
        # Real variables
        sim_end=merge(sim_end,(X_Laspeyres = sum(sim_end.measure.*p_o_0.ξ.*p_o_0.pf.*p_o_end.yf),
        D_Laspeyres = sum(sim_end.measure.* p_o_0.pd.*p_o_end.yd),
        Sales_Laspeyres = sum(sim_end.measure.*(p_o_0.pd.*p_o_end.yd + p_o_0.ξ.*p_o_0.pf.*p_o_end.yf)),
        GDP_Laspeyres = sum(sim_end.measure.*(p_o_0.pd.*p_o_end.yd + p_o_0.ξ*p_o_0.pf.*p_o_end.yf -p_o_end.m*p_o_0.Pk))))

 # STEP 3: guess sequence of aggregate prices for N period, with N sufficiently large [this was in a separate script in matlab, KLS3_trans_obj, you could that here but there isn't much gain and it's useful to mantain a low number of scripts]

 # Store results from initial and final steady states

rt=Vector{Any}(undef,p.N) #this creates a vector that accept every type in its entries (hence the "any")

rt[1] = merge(p_o_0,(measure=sim_0.measure,))
rt[p.N] = merge(p_o_end,(measure=sim_end.measure,))
 #r_0=nothing #[Julia's equivalent to clear]
 #r_end=nothing

 # Guess sequence of aggregate prices for N period, with N sufficiently large
 #ϕht = zeros(p.N,1)
 p_o_t=(Pkt = zeros(p.N,1),
 wt = zeros(p.N,1),
 ξt = zeros(p.N,1),
 Ykt = zeros(p.N,1),
 Yct = zeros(p.N,1))

    if s.transition_AD == 1 # If AutoDiff is enabled, tariffsincomet must be pre-allocated in a way that accepts Real numbers (and not just Float64 as default). This is because when the solver uses AD, it uses a different number type (DualNumbers).
    p_o_t=merge(p_o_t,(tariffsincomet = zeros(Real,p.N,1),))
    else
    p_o_t=merge(p_o_t,(tariffsincomet = zeros(p.N,1),))
    end

    if s.tariffsincome == 1
     # Guess
     #p_o_t.ϕht[1] = p_o_0.ϕ_h
     p_o_t.wt[1] = p_o_0.w
     p_o_t.ξt[1] = p_o_0.ξ
     p_o_t.Pkt[1] = p_o_0.Pk
     p_o_t.Yct[1] = p_o_0.Yc
     p_o_t.Ykt[1] = p_o_0.Yk

        if s.PE == 0
         #p_o_t.ϕht[p.N] = p_o_end.ϕ_h
         p_o_t.wt[p.N] = p_o_end.w
         p_o_t.ξt[p.N] = p_o_end.ξ
         p_o_t.Pkt[p.N] = p_o_end.Pk
         p_o_t.Yct[p.N] = p_o_end.Yc
         p_o_t.Ykt[p.N] = p_o_end.Yk

        elseif s.PE ==1

         #p_o_t.ϕht[p.N]=p_o_0.ϕ_h
         p_o_t.wt[p.N] = p_o_0.w
         p_o_t.ξt[p.N] = p_o_0.ξ
         p_o_t.Pkt[p.N]=p_o_end.Pk #This was in the MATLAB code, is it a typo and should it be p_o_0.Pk as the rest of the prices?
         p_o_t.Yct[p.N] = p_o_0.Yc
         p_o_t.Ykt[p.N] = p_o_0.Yk

        end

     const period2_change = 0.5
     #I labeled every init/expo differently because as I am using "const" (see lines above for an explanation), using the same name for all would be problematic (because it would make them not constant after all)
     const init1 =p_o_t.wt[1]+ (p_o_t.wt[p.N]-p_o_t.wt[1])*period2_change
     const expo1 =log.(p_o_t.wt[p.N]./init1)./log.(p.N)   #p_o_t.wt[1]*(p.N)^expo1=p_o_t.wt[p.N] =>expo1=log.(p_o_t.wt[p.N]/p_o_t.wt[1])/log.(p.N)
     p_o_t.wt[2:p.N-1].= @. init1*(2:p.N-1).^expo1

     const init2 =p_o_t.Yct[1]+ (p_o_t.Yct[p.N]-p_o_t.Yct[1])*period2_change
     const expo2 =log.(p_o_t.Yct[p.N]./init2)./log.(p.N) #p_o_t.wt[1]*(p.N)^expo=p_o_t.wt[p.N] =>expo2=log.(p_o_t.wt[p.N]/p_o_t.wt[1])/log.(p.N)
     p_o_t.Yct[2:p.N-1].= @. init2*(2:p.N-1).^expo2

     const init3 =p_o_t.Ykt[1]+ (p_o_t.Ykt[p.N]-p_o_t.Ykt[1])*period2_change
     const expo3 =log.(p_o_t.Ykt[p.N]./init3)./log.(p.N)   #p_o_t.wt[1]*(p.N).^expo3=p_o_t.wt[p.N] =>expo3=log.(p_o_t.wt[p.N]./p_o_t.wt[1])/log.(p.N)
     p_o_t.Ykt[2:p.N-1].= @. init3*(2:p.N-1).^expo3

     p_o_t.Pkt[2:p.N-1].= p_o_end.Pk
     p_o_t.ξt[2:p.N-1].=p_o_end.ξ

     Guess = [p_o_t.Yct[2:p.N-1] p_o_t.wt[2:p.N-1] p_o_t.ξt[2:p.N-1] p_o_t.Pkt[2:p.N-1] p_o_t.Ykt[2:p.N-1]]

        if s.load_prices == 1
        @load "KLS3_Prices.jld2" Prices
         Guess = [Prices[1,2:p.N-1] Prices[2,2:p.N-1] Prices[3,2:p.N-1] Prices[4,2:p.N-1] Prices[5,2:p.N-1]]
        end
    end


## STEP 4: Solving for the GE prices and quantities

optim_trans(F,x)=KLS3_transition_vec2!(F,x,p,p_o_t,rt) #"instantiating" the objective functions with the previously specified values of the parameters/pre-allocations (in order to use it as an argument for the solvers below)

@time begin
        if s.transition_GE == 1 && s.solver_LM==0
            if s.transition_AD == 1
                #using NLsolve (Trust Region Dogleg Algorithm)
                results_trans = nlsolve(optim_trans,log.(Guess), autodiff = :forward, method=s.method_trans, factor=1, autoscale=false, xtol=s.xtol_trans, ftol=s.ftol_trans,iterations=s.MaxIter_trans, show_trace=s.show_trace_trans,store_trace=true)
            else
                results_trans = nlsolve(optim_trans,log.(Guess), method=s.method_trans, factor=1, autoscale=true, xtol=s.xtol_trans, ftol=s.ftol_trans,iterations=s.MaxIter_trans, show_trace=s.show_trace_trans,store_trace=true)
            end
            Prices_sol=results_trans.zero
            mc, p_o_t, sim_fun, rt = KLS3_transition_vec2(Prices_sol,p,p_o_t,rt)
        elseif s.transition_GE == 1 && s.solver_LM==1
            #using LeastSquaresOptim #(L-M algorithm)
            # #[This solver allows the use of the Levenberg-Marquardt algorithm, which seems to be better at computing the transition (it throws out almost the same solution as matlab's, while NLsolve throws a much noiser one; this changes if the guess is closer to the solution, though). Unfortunately, it takes around 10x the time NLsolve takes] IMPORTANT: If using this solver, uncomment first lines of function KLS3_transition_vec2!

            if s.transition_AD == 1
                solver_LSO=optimize!(LeastSquaresProblem(x = vec(log.(Guess)), f! = optim_trans, output_length = length(vec(log.(Guess))), autodiff = :forward), x_tol=s.xtol_trans, f_tol=s.ftol_trans,iterations=s.MaxIter_trans,LevenbergMarquardt())
            else
                solver_LSO=optimize!(LeastSquaresProblem(x = vec(log.(Guess)), f! = optim_trans, output_length = length(vec(log.(Guess))), autodiff = :forward), x_tol=s.xtol_trans, f_tol=s.ftol_trans,iterations=s.MaxIter_trans,LevenbergMarquardt())
            end

            Prices_sol=hcat(solver_LSO.minimizer[1:p.N-2],solver_LSO.minimizer[p.N-1:2*p.N-4],solver_LSO.minimizer[2*p.N-3:3*p.N-6],solver_LSO.minimizer[3*p.N-5:4*p.N-8],solver_LSO.minimizer[4*p.N-7:5*p.N-10])
            mc_lso, p_o_t_lso, sim_fun_lso, rt_lso = KLS3_transition_vec2(Prices_sol,p,p_o_t,rt)
        else
            Prices_sol = log.(Guess)
            mc, p_o_t, sim_fun, rt = KLS3_transition_vec2(Prices_sol,p,p_o_t,rt)
        end
    say("Transition Finished - you'd better come and take a look....")
end

end

 ## Welfare analysis #INCOMPLETE [3D FIGURES PENDING]
 if s.welfare==1
     include("welfare.jl")
 end

 ## End

 if s.save_prices==1
     if s.tariffsincome == 0
         Prices[1,:] = [m.ϕht[1] exp.(Prices_sol[1:p.N-2])' m.ϕht[p.N]]
         Prices[2,:] = [m.wt[1] exp.(Prices_sol[p.N-1:2*p.N-4])' m.wt[p.N]]
         Prices[3,:] = [m.ξt[1] exp.(Prices_sol[2*p.N-3:3*p.N-6])' m.ξt[p.N]]
         Prices[4,:] = [m.Pkt[1] exp.(Prices_sol[3*p.N-5:4*p.N-8])' m.Pkt[p.N]]
     elseif s.tariffsincome == 1
         Prices[1,:] = [m.Yct[1] exp.(Prices_sol[1:p.N-2])' m.Yct[p.N]]
         Prices[2,:] = [m.wt[1] exp.(Prices_sol[p.N-1:2*p.N-4])' m.wt[p.N]]
         Prices[3,:] = [m.ξt[1] exp.(Prices_sol[2*p.N-3:3*p.N-6])' m.ξt[p.N]]
         Prices[4,:] = [m.Pkt[1] exp.(Prices_sol[3*p.N-5:4*p.N-8])' m.Pkt[p.N]]
         Prices[5,:] = [m.Ykt[1] exp.(Prices_sol[4*p.N-7:end])' m.Ykt[p.N]]
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
