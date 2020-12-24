



##################### Directory settings #####################
    results=@with_kw (dir = pwd(),
    dir_results = string(results,"/","results", "/"))
##################### Tariff Income #####################
    t_inc=@with_kw (tariffsincome = 1, # = 0 if don't give tariffs income back
                             # = 1 give tariffs income back (we need additional
                             # guess)
    tariffsincome_PE_flag=1,
    tariffsincome_PE = 0)

##################### Steady State in GE? #####################
    ge=@with_kw (GE =1, # =0 if partial equilibrium (inital SS)
                  # =1 if general equilibrium (initial SS)
        GE_end = 1) # =0 if partial equilibrium (final SS)
                      # =1 if general equilibrium (final SS)

##################### Additional Options #####################

#Denomination of fixed costs (includes sunk costs)
# =0 if fixed costs denominated in units of labor
# =1 if fixed costs denominated in units of the final good
    ad_opt=@with_kw (fcost_fgoods = 0,

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
    load_prices = 0) # load prices


##################### Baseline Parameters #####################

p_dft = @with_kw (θ = 0.20541136, #Collateral constraint
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
log_z_μ = log(z_mean)-(log_z_σ^2)*(1/((1-log_z_ρ^2)))*(1/2)) #Normalize average productivity to 1, log_z_mu is the mean of log_z


##################### Transition Settings #####################

m=parameter_defaults()
N=60 #length of transition

#Shock to collateral constraint
θ_old = m.θ;
θ_new = m.θ;
θ_v = ones(N,1);
θ_v[1:2] = θ_old;
θ_v[3:end] = θ_new;

#Real interest rate shock
r_old = m.r;
r_new = m.r;
rv = zeros(N,1);
rv[1:2] = r_old;
rv[3:end] = r_new;

# Beta shock
β_old = m.β;
β_new = m.β;
β_v = zeros(N,1);
β_v[1] = β_old;
β_v[2:end] = β_new;

# delta shock
δ_old = m.δ;
δ_new = m.δ;
δ_v = zeros(N,1);
δ_v[1] = δ_old;
δ_v[2:end] = δ_new;

# Foreign CPI shock
Pf_old = m.Pf;
Pf_new = m.Pf;
Pfv = zeros(N,1);
Pfv[1] = Pf_old;
Pfv[2:end] = Pf_new;

#Shock to pm
# pm_old = m.Pm;
# pm_new = m.Pm;
# pm_v = zeros(N,1);
# pm_v[1] = pm_old;
# pm_v[2:end] = pm_new;

# Foreign output shock
# Negative shock -> devaluation
# Positive shock -> apreciation
Yf_old = m.Yf;
Yf_new = m.Yf;
Yfv = zeros(N,1);
Yfv[1] = Yf_old;
Yfv[2:end] = Yf_new;

#Shock to iceberg costs
τ_old = m.τ;
τ_new = 1*m.τ;
τ_v = zeros(N,1);
τ_v[1:2] = τ_old;
τ_v[3:end] = τ_new;

#Shocks to tariffs
τ_m_c_old = m.τ_m_c;
τ_m_c_new = 0.375*m.τ_m_c;
τ_m_c_v = zeros(N,1);
τ_m_c_v[1] = τ_m_c_old;
τ_m_c_v[2:end] = τ_m_c_new;

τ_m_k_old = m.τ_m_k;
τ_m_k_new = 1*m.τ_m_k;
τ_m_k_v = zeros(N,1);
τ_m_k_v[1] = τ_m_k_old;
τ_m_k_v[2:end] = τ_m_k_new;


τ_x_old = m.τ_x;
τ_x_new = 1*m.τ_x;
τ_x_v = zeros(N,1);
τ_x_v[1] = τ_x_old;
τ_x_v[2:end] = τ_x_new;

t_s = @with_kw (N=N,
θ_v=θ_v,
rv=rv,
β_v=β_v,
δ_v=δ_v,
Pfv=Pfv,
Yfv=Yfv,
τ_v=τ_v,
τ_m_c_v=τ_m_c_v,
τ_m_k_v=τ_m_k_v,
τ_x_v=τ_x_v)

##################### Solution Options #####################

s_op=@with_kw (
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
     MaxIter = 15)

     op_GE =@with_kw (method=:trustregion,xtol=1E-7,ftol=1E-8,show_trace=true)

     # This has to be translated to NLsolve options [PENDING]
     options_trans = optimoptions('fsolve','Display','iter','StepTolerance',1e-7,'FunctionTolerance',1e-7,'MaxFunctionEvaluations',s.MaxFunEvalsTrans,'MaxIterations',s.MaxIter,...
     'Algorithm','levenberg-marquardt','InitDamping',0.01,'FiniteDifferenceStepSize',0.001,'FiniteDifferenceType','central'...
     ,'ScaleProblem','jacobian','UseParallel',true)

     ## [This is not translated, copy/paste from Matlab] Simulation by generating sequence of shocks (if s.flag_simulate==1)

     # if ad_opt.flag_simulate == 1
     #     #Simulations
     #     sim_op=@withkw(T = 40, #50;   # Number of periods
     #     Nfirms = 1000000, #250000; #1000000; #400000;  # Number of firms
     #     burn = 1, # Number of time periods burnt = s.burn*s.T
     #     seed = 88,
    #      setGlobalStream(RandStream('mt19937ar','seed',s.seed)); # Sets seed for randomizations)


    ##################### Setup asset grids, productivity grids, annuity market transfers #####################

        [log_z_grid, z_P, z_π] = tauchen(s_op.z_grid_size,4,p_dft.log_z_ρ,p_dft.log_z_σ,p_dft.log_z_μ,s_op.z_grid_power)
        z_grid = exp(log_z_grid');
        z_grid_original = exp(log_z_grid');

        grids=@with_kw (log_z_grid=log_z_grid, #Productivity
        z_P=z_P,
        z_π=z_π,
        z_grid=exp(log_z_grid'),
        z_grid_original=exp(log_z_grid'),
        z_min=min(z_grid), #Productivity process statistics
        z_max = max(z_grid),
        z_mean = sum(z_grid.*z_π),
        z = ones(s_op.a_grid_size,1)*z_grid',
        a_grid = LinRange(0,1,s_op.a_grid_size), #Asset grid
        a_grid = a_grid.^s_op.a_grid_power,
        a_grid = a_grid.*(s_op.a_grid_ub-s_op.a_grid_lb) + s_op.a_grid_lb,
        a_grid_vec = repeat(a_grid,1,length(a_grid)))
