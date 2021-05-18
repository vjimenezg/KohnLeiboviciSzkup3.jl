
### Tauchen Function, Parameters Setup and Settings ###

### Tauchen ###

function tauchen1(n::Integer,ρ::Real, σ_e::Real, μ=zero(typeof(ρ)),c::Real=4,z_power::Real=1)

#Standard deviation of stationaty process y (and discrete process z)
σ_z=σ_e/sqrt(1-ρ^2)

#Range
r=(2*σ_z)*c

#Solve for values of discrete state space
z_1 = μ - c*σ_z
z_end = μ + c*σ_z
z = (LinRange(0,1,n).^(z_power)).*(z_end-z_1) .+ z_1
z = z'

#Construct midpoints m(i), i=1...n-1
m = (z[2:n]+z[1:n-1])/2

P = zeros(n,n)
d=Distributions.Normal()
for i=1:n
P[i,1] = Distributions.cdf(d,(m[1]-(1-ρ)*μ-ρ*z[i])/σ_e)
   for j=2:n-1
        P[i,j] = Distributions.cdf(d,(m[j]-(1-ρ)*μ-ρ*z[i])/σ_e)-Distributions.cdf(d,(m[j-1]-(1-ρ)*μ-ρ*z[i])/σ_e)
   end
P[i,n] = 1-Distributions.cdf(d,(m[n-1]-(1-ρ)*μ-ρ*z[i])/σ_e)
end

mc=QuantEcon.MarkovChain(P)
π_1=hcat(QuantEcon.stationary_distributions(mc)...)'

return z,P,π_1
end

### Parameters Setup ###

#What @with_kw (Parameters) does is to create a struct (new type) that contains named elements with default values. This allows a better performance (as opposed to using a regular tuple). This object can be instantiated through parameters_default() with the default values or through parameters_default(name=xx) if you want to set up the parameter name to a value xx different from the default one.
parameters_default = @with_kw (θ = 0.20541136, #Collateral constraint
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
ω_h_c = 1.0,
ω_m_c = 0.20904358,
ω_h_k = 1.0,
ω_m_k = 0.28722687,
Yf = 10.0, #/((1+m.τ_x)^(-m.σ)) #Output in the rest of the world
Pf = 1.0, #Price index in the rest of the world
Pm_c = 1.0, #/(1+m.τ_m); #Price index of goods imported from the rest of the world
Pm_k = 1.0,
ω_m = 1.0, #Measure of varieties produced by the rest of the world
r = 0.06,
tariffsincome = 0.0,
log_z_ρ = 0.9,
log_z_μ = log(z_mean)-(log_z_σ^2)*(1/(1.0-log_z_ρ^2))*(1/2), #Normalize average productivity to 1, log_z_μ is the mean of log_z
N=20,
θ_old = θ,
θ_new = θ,
θ_v = hcat(θ_old*ones(1,2),θ_new*ones(1,N-2)),
r_old = r,
r_new = r,
rv = hcat(r_old*ones(1,2),r_new*ones(1,N-2)),
β_old = β,
β_new = β,
β_v = hcat(β_old*ones(1,1),β_new*ones(1,N-1)),
δ_old = δ,
δ_new = δ,
δ_v = hcat(δ_old*ones(1,1),δ_new*ones(1,N-1)),
Pf_old = Pf,
Pf_new = Pf,
Pfv = hcat(Pf_old*ones(1,1),Pf_new*ones(1,N-1)),
Yf_old = Yf,
Yf_new = Yf,
Yfv = hcat(Yf_old*ones(1,1),Yf_new*ones(1,N-1)),
τ_old = τ,
τ_new = 1*τ,
τ_v = hcat(τ_old*ones(1,2),τ_new*ones(1,N-2)),
τ_m_c_old = τ_m_c,
τ_m_c_new = 1*τ_m_c,
τ_m_c_v = hcat(τ_m_c_old*ones(1,1),τ_m_c_new*ones(1,N-1)),
τ_m_k_old = τ_m_k,
τ_m_k_new = 0.12, #0.27
τ_m_k_v = hcat(τ_m_k_old*ones(1,1),τ_m_k_new*ones(1,N-1)),
τ_x_old = τ_x,
τ_x_new = 1*τ_x,
τ_x_v = hcat(τ_x_old*ones(1,1),τ_x_new*ones(1,N-1)),
z_grid_size = 100, # 100; #75; #250;# Productivity grid size
z_grid_power =1/2, # 1/2; #1/2; #1; # Curvature parameter to control
                                    # distance across grid points

#Assets
a_grid_size = 100,# 100; #150; #250; # Asset grid size
a_grid_power = 2, #2; #3 # Curvature parameter to control distance across
                         # grid points -- for a>=0 (for a<0, grid spacing is linear)

a_grid_ub = 5, #200; #200; #1e-3 #500;# Upper bound on asset grid
                                      # CHECK WHETHER IT BINDS!
a_grid_lb = 1e-3,
disc =tauchen1(z_grid_size,log_z_ρ,log_z_σ,log_z_μ,4,z_grid_power),
log_z_grid=disc[1],
z_P=disc[2],
z_π=disc[3]',
z_grid = exp.(log_z_grid'),
z_grid_original = exp.(log_z_grid'),
z_min=minimum(z_grid,dims=1),
z_max = maximum(z_grid,dims=1),
z_mean1 = sum(z_grid.*z_π),
z = ones(a_grid_size,1)*z_grid',

#Asset grid
a_grid_1 = LinRange(0,1,a_grid_size),
a_grid_2 = a_grid_1.^a_grid_power,
a_grid = (a_grid_2.*(a_grid_ub-a_grid_lb) .+ a_grid_lb)',
a_grid_vec = repeat(a_grid,length(a_grid),1))

### Settings ###
#same as parameters_default, a struct with named variables inside.

settings_default = @with_kw (dir=pwd(),
dir_results = string(dir,"/","results", "/"),
tariffsincome = 1, # = 0 if don't give tariffs income back # = 1 give tariffs income back (we need additional guess)
tariffsincome_PE_flag=1,
solver_LM=0, # = 1 if LeastSquaresOptim with Levenberg-Marquardt Algorithm, = 0 if NLsolve with Trust Region Algorithm
tariffsincome_PE = 0,
GE =1, # =0 if partial equilibrium (inital SS)
       # =1 if general equilibrium (initial SS)
GE_end = 1, # =0 if partial equilibrium (final SS)
            # =1 if general equilibrium (final SS))
fcost_fgoods = 0,
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
transition_AD=1, #additional option, to use or not to use automatic differentiation (see https://julia.quantecon.org/more_julia/optimization_solver_packages.html#Introduction-to-Differentiable-Programming for an explanation). This shouldn't be turned off unless there are issues with the code, autodifferentiation makes the code (weakly) faster and more accurate.
PE = 0, # transition and final SS using initial SS prices
welfare = 1, # Experiments
save_workspace = 1, # Saving options
save_prices = 1,
load_prices = 0,
e_speedup = 1,

#Convergence parameters
eps = 1e-8,
eps_sm = 1e-10, #Stationary measure

ac_iter = 10, #8; #Accelerator iterations

#Optimizer options

 method_GE=:trust_region, #Method used in the solver, I tried the other two available in NLsolve and this is the best in terms of accuracy (though anderson can sometimes be faster)
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
 MaxIter_trans = 15)
