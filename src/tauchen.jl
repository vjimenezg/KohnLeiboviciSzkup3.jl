function tauchen(n::Integer,ρ::Real, σ_e::Real, μ=zero(typeof(ρ)),c::Real=4,z_power::Real=1)

σ_z=σ_e/sqrt(1-ρ^2)
r=(2*σ_z)*c
z_1 = μ - c*σ_z
z_end = μ + c*σ_z
z = (LinRange(0,1,n).^(z_power)).*(z_end-z_1) .+ z_1
z = z'

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
π_1=QuantEcon.stationary_distributions(mc)

return z,P,π_1
end
