

##################### Directory settings #####################
    results = (dir = pwd(),)
    results=merge(results,(dir_results = string(results.dir,"/","results", "/"),))
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
 #Welfare Graphs
 # =0 if graphs are to be done directly in Julia (not available yet)
 # =1 if welfare data is to be exported to matlab
    graphs_matlab = 1,

# Simulation
    flag_simulate = 0, # flag_simulate = 0, simulates by iterating on measure
                       # of firms, flag_simulate = 1, simulates by generating
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

m = (θ = 0.20541136, #Collateral constraint
    F_base = 0.46563816, #Fixed export cost
    log_z_σ = 0.15652825, #Standard dev. of log-normal productivity distribution
α_m=0.5,
β=0.83500441,
σ=4.0,
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
Yf = 10, #/((1+m.τ_x)^(-m.σ)) #Output in the rest of the world
Pf = 1.0, #Price index in the rest of the world
Pm_c = 1.0, #/(1+m.τ_m); #Price index of goods imported from the rest of the world
Pm_k = 1.0,
ω_m = 1.0, #Measure of varieties produced by the rest of the world
r = 0.06,
tariffsincome = 0.0,
log_z_ρ = 0.9) #13.69641182/14.69641182, #Persistence
m=merge(m,(log_z_μ = log(m.z_mean)-(m.log_z_σ^2)*(1/(1-m.log_z_ρ^2))*(1/2),))
#Normalize average productivity to 1, log_z_μ is the mean of log_z


##################### Transition Settings #####################
s=merge(s,(N=20,  #length of transition
transition_AD=0)) #Autodifferentiation on/off

#Shock to collateral constraint
θ_old = m.θ
θ_new = m.θ
m=merge(m,(θ_v = hcat(θ_old*ones(1,2),θ_new*ones(1,s.N-2)),))

#Real interest rate shock
r_old = m.r
r_new = m.r
m=merge(m,(rv = hcat(r_old*ones(1,2),r_new*ones(1,s.N-2)),))


# β shock
β_old = m.β
β_new = m.β
m=merge(m,(β_v = hcat(β_old*ones(1,1),β_new*ones(1,s.N-1)),))


# δ shock
δ_old = m.δ
δ_new = m.δ
m=merge(m,(δ_v = hcat(δ_old*ones(1,1),δ_new*ones(1,s.N-1)),))

# Foreign CPI shock
Pf_old = m.Pf
Pf_new = m.Pf
m=merge(m,(Pfv = hcat(Pf_old*ones(1,1),Pf_new*ones(1,s.N-1)),))


#Shock to pm
# pm_old = m.Pm
# pm_new = m.Pm
# m=merge(m,(pm_v = hcat(pm_old*ones(1,1),pm_new*ones(1,s.N-1)),))

# Foreign output shock
# Negative shock -> devaluation
# Positive shock -> apreciation
Yf_old = m.Yf
Yf_new = m.Yf
m=merge(m,(Yfv = hcat(Yf_old*ones(1,1),Yf_new*ones(1,s.N-1)),))


#Shock to iceberg costs
τ_old = m.τ
τ_new = 1*m.τ
m=merge(m,(τ_v = hcat(τ_old*ones(1,2),τ_new*ones(1,s.N-2)),))


#Shocks to tariffs
τ_m_c_old = m.τ_m_c
τ_m_c_new = 1*m.τ_m_c
m=merge(m,(τ_m_c_v = hcat(τ_m_c_old*ones(1,1),τ_m_c_new*ones(1,s.N-1)),))


τ_m_k_old = m.τ_m_k
τ_m_k_new = 0.12 #0.27
m=merge(m,(τ_m_k_v = hcat(τ_m_k_old*ones(1,1),τ_m_k_new*ones(1,s.N-1)),))


τ_x_old = m.τ_x
τ_x_new = 1*m.τ_x
m=merge(m,(τ_x_v = hcat(τ_x_old*ones(1,1),τ_x_new*ones(1,s.N-1)),))


##################### Solution Options #####################

s=merge(s,(
#Productivity
    z_grid_size = 100, # 100; #75; #250;# Productivity grid size
    z_grid_power =1/2, # 1/2; #1/2; #1; # Curvature parameter to control
                                        # distance across grid points

#Assets
    a_grid_size = 100,# 100; #150; #250; # Asset grid size
    a_grid_power = 2, #2; #3 # Curvature parameter to control distance across
                             # grid points -- for a>=0 (for a<0, grid spacing is linear)

    a_grid_ub = 5, #200; #200; #1e-3 #500;# Upper bound on asset grid
                                          # CHECK WHETHER IT BINDS!
    a_grid_lb = 1e-3, #1e-8;

#Export entry solution algorithm
#Constrain export policy to be monotonic for assets (given productivity)
    e_speedup = 1,

#Convergence parameters
    eps = 1e-8,
    eps_sm = 1e-10, #Stationary measure

    ac_iter = 10, #8; #Accelerator iterations

#Optimizer options

     method_GE=:trust_region,
     xtol_GE=1E-7,
     ftol_GE=1E-8,
     show_trace_GE=false,
     #MaxFunEvals_GE = 8_000,
     #MaxIter_GE = 15,

     method_trans=:trust_region,
     xtol_trans=1E-7,
     ftol_trans=1E-7,
     show_trace_trans=true,
     MaxFunEvals_trans = 8_000,
     MaxIter_trans = 15))



    ##################### Setup asset grids, productivity grids, annuity market transfers #####################

        #Productivity

        log_z_grid,z_P,z_π =tauchen(s.z_grid_size,m.log_z_ρ,m.log_z_σ,m.log_z_μ,4,s.z_grid_power)

        r = (log_z_grid=log_z_grid, z_P=z_P, z_π=z_π')
        r=merge(r,(z_grid = exp.(r.log_z_grid'),
        z_grid_original = exp.(r.log_z_grid')))

        #Productivity process statistics
        r=merge(r,(z_min=minimum(r.z_grid,dims=1),
        z_max = maximum(r.z_grid,dims=1),
        z_mean = sum(r.z_grid.*r.z_π),
        z = ones(s.a_grid_size,1)*r.z_grid',

        #Asset grid
        a_grid_1 = LinRange(0,1,s.a_grid_size))) #Asset grid
        r=merge(r,(a_grid_2 = r.a_grid_1.^s.a_grid_power,))
        r=merge(r,(a_grid = (r.a_grid_2.*(s.a_grid_ub-s.a_grid_lb) .+ s.a_grid_lb)',))
        r=merge(r,(a_grid_vec = repeat(r.a_grid,length(r.a_grid),1),))
