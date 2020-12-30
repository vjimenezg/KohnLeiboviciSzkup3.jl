



##################### Directory settings #####################
    results=(dir = pwd(),
    dir_results = string(dir,"/","results", "/"))
##################### Tariff Income #####################
    s=(tariffsincome = 1, # = 0 if don't give tariffs income back
                             # = 1 give tariffs income back (we need additional
                             # guess)
    tariffsincome_PE_flag=1,
    tariffsincome_PE = 0)

##################### Steady State in GE? #####################
    s=merge(s,(GE =1, # =0 if partial equilibrium (inital SS)
                  # =1 if general equilibrium (initial SS)
        GE_end = 1)) # =0 if partial equilibrium (final SS)
                      # =1 if general equilibrium (final SS)

##################### Additional Options #####################

#Denomination of fixed costs (includes sunk costs)
# =0 if fixed costs denominated in units of labor
# =1 if fixed costs denominated in units of the final good
    s=merge(s,(fcost_fgoods = 0,

# Simulation
    flag_simulate = 0, # flag_simulate = 0, simulates by iterating on measure of firms,
                         # flag_simulate = 1, simulates by generating
                         # sequence of random shocks
    display = 1, #display=1, displays GE results on each iteration, 0 does not.
    extra_results=0,
    transition = 1, #Transition
    transition_GE = 1,
    PE = 0, # transition and finanl SS using initial SS prices
    welfare = 1, # Experiments
    save_workspace = 1, # Saving options
    save_prices = 1,
    load_prices = 0)) # load prices


##################### Baseline Parameters #####################

m= (θ = 0.20541136, #Collateral constraint
    F_base = 0.46563816, #Fixed export cost
    log_z_σ = 0.15652825, #Standard deviation of log-normal productivity distribution
α_m=0.5,
β=0.83500441,
σ = 4.0,
γ=2.0,
α=0.6,
δ=0.1,
z_mean=1.0,
τ=1.0,
τ_m_c = 0.32,
τ_m_k = 0.32,
τ_x = 0.32,
ω_h_c = 1,
ω_m_c = 0.20904358,
ω_h_k = 1,
ω_m_k = 0.28722687,
Yf = 10, #/((1+m.tau_x)^(-m.sigma)); #Output in the rest of the world
Pf = 1.0, #Price index in the rest of the world
Pm_c = 1.0, #/(1+m.tau_m); #Price index of goods imported from the rest of the world
Pm_k = 1.0,
ω_m = 1.0, #Measure of varieties produced by the rest of the world # CHECK)
r = 0.06,
tariffsincome = 0.0,
log_z_ρ = 0.9, #13.69641182/14.69641182, #Persistence
log_z_μ = log(m.z_mean)-(m.log_z_σ^2)*(1/((1-m.log_z_ρ^2)))*(1/2)) #Normalize average productivity to 1, log_z_mu is the mean of log_z


##################### Transition Settings #####################

s=merge(s,(N=60,)) #length of transition

#Shock to collateral constraint
θ_old = m.θ;
θ_new = m.θ;
m=merge(m,(θ_v = vcat(θ_old*ones(2,1),θ_new*ones(N-2,1)),))

#Real interest rate shock
r_old = m.r;
r_new = m.r;
m=merge(m,(rv = vcat(r_old*ones(2,1),r_new*ones(N-2,1)),))


# Beta shock
β_old = m.β;
β_new = m.β;
m=merge(m,(β_v = vcat(β_old,β_new*ones(N-1,1)),))


# delta shock
δ_old = m.δ;
δ_new = m.δ;
m=merge(m,(δ_v = vcat(δ_old,δ_new*ones(N-1,1)),))

# Foreign CPI shock
Pf_old = m.Pf;
Pf_new = m.Pf;
m=merge(m,(Pfv = vcat(Pf_old,Pf_new*ones(N-1,1)),))


#Shock to pm
# pm_old = m.Pm;
# pm_new = m.Pm;
# m=merge(m,(pm_v = vcat(pm_old,pm_new*ones(N-1,1)),))

# Foreign output shock
# Negative shock -> devaluation
# Positive shock -> apreciation
Yf_old = m.Yf;
Yf_new = m.Yf;
m=merge(m,(Yfv = vcat(Yf_old,Yf_new*ones(N-1,1)),))


#Shock to iceberg costs
τ_old = m.τ;
τ_new = 1*m.τ;
m=merge(m,(τ_v = vcat(τ_old*ones(2,1),τ_new*ones(N-2,1)),))


#Shocks to tariffs
τ_m_c_old = m.τ_m_c;
τ_m_c_new = 0.375*m.τ_m_c;
m=merge(m,(τ_m_c_v = vcat(τ_m_c_old,τ_m_c_new*ones(N-1,1)),))


τ_m_k_old = m.τ_m_k;
τ_m_k_new = 1*m.τ_m_k;
m=merge(m,(τ_m_k_v = vcat(τ_m_k_old,τ_m_k_new*ones(N-1,1)),))


τ_x_old = m.τ_x;
τ_x_new = 1*m.τ_x;
m=merge(m,(τ_x_v = vcat(τ_x_old,τ_x_new*ones(N-1,1)),))


##################### Solution Options #####################

s=merge(s,(
#Productivity
    z_grid_size = 200, # #100; #75; #250; #Productivity grid size
    z_grid_power =1/2, # 1/2; #1/2; #1; #Curvature parameter to control distance across grid points

#Assets
    a_grid_size = 200,# 100; #150; #250; #Asset grid size
    a_grid_power = 2, #2; #3 #Curvature parameter to control distance across grid points -- for a>=0 (for a<0, grid spacing is linear)

    a_grid_ub = 5, #200; #200; #1e-3 #500;#Upper bound on asset grid #CHECK WHETHER IT BINDS!
    a_grid_lb = 1e-3, #1e-8;

#Export entry solution algorithm
#Constrain export policy to be monotonic for assets (given productivity)
    e_speedup = 1,

#Convergence parameters
    eps = 1e-8,
    eps_sm = 1e-10, #Stationary measure

    ac_iter = 10, #8; #Accelerator iterations

#Optimizer options
     MaxFunEvalsTrans = 8000,
     MaxIter = 15,

     method_GE=:trustregion,
     xtol_GE=1E-7,
     ftol_GE=1E-8,
     show_trace_GE=true))

     # This has to be translated to NLsolve options [PENDING]

    # options_trans = optimoptions('fsolve','Display','iter','StepTolerance',1e-7,'FunctionTolerance',1e-7,'MaxFunctionEvaluations',s.MaxFunEvalsTrans,'MaxIterations',s.MaxIter,...
    # 'Algorithm','levenberg-marquardt','InitDamping',0.01,'FiniteDifferenceStepSize',0.001,'FiniteDifferenceType','central'...
    # ,'ScaleProblem','jacobian','UseParallel',true)

     ## [This is not translated, copy/paste from Matlab, originally commented] Simulation by generating sequence of shocks (if s.flag_simulate==1)

     # if ad_opt.flag_simulate == 1
     #     #Simulations
     #     sim_op=@withkw(T = 40, #50;   # Number of periods
     #     Nfirms = 1000000, #250000; #1000000; #400000;  # Number of firms
     #     burn = 1, # Number of time periods burnt = s.burn*s.T
     #     seed = 88,
    #      setGlobalStream(RandStream('mt19937ar','seed',s.seed)); # Sets seed for randomizations)


    ##################### Setup asset grids, productivity grids, annuity market transfers #####################

        #Productivity

        log_z_grid,z_P,z_π =tauchen(s.z_grid_size,4,m.log_z_ρ,m.log_z_σ,m.log_z_μ,s.z_grid_power)

        r=(log_z_grid=log_z_grid, z_P=z_P, z_π=z_π,
        z_grid = exp(r.log_z_grid'),
        z_grid_original = exp(r.log_z_grid'),

        #Productivity process statistics
        z_min=min(r.z_grid),
        z_max = max(r.z_grid),
        z_mean = sum(r.z_grid.*r.z_π),
        z = ones(s.a_grid_size,1)*r.z_grid',

        #Asset grid
        a_grid = LinRange(0,1,s.a_grid_size), #Asset grid
        a_grid = r.a_grid.^s.a_grid_power,
        a_grid = r.a_grid.*(s.a_grid_ub-s.a_grid_lb) + s.a_grid_lb,
        a_grid_vec = repeat(r.a_grid,1,length(r.a_grid)))
