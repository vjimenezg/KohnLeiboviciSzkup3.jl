#Baseline parameters

parameter_defaults = @with_kw (θ = 0.20541136, #Collateral constraint
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
