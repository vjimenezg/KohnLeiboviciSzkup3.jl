function KLS3_GE_par(x,m,s,r)

#### Aggregate prices and quantities ####

#Guessed prices

w = exp(x[1])
m.Φ_h = exp(x[2])
m.ξ = exp(x[3])
m.Pk = exp(x[4])

m.Pk_lag = m.Pk


### Solution

# Tariff income (initialized to 0)
if s.tariffsincome == 1
    m.Yk = exp(x[5])
    m.Yc = exp(x[2])
    m.Φ_h = (m.ω_h_c^m.σ)*m.Yc + ((m.ω_h_k*m.Pk)^m.σ)*m.Yk
    %Yc =  (m.Φ_h - ( m.ω_m_k*m.Pk)^m.σ*m.Yk)/( m.ω_m_c^m.σ)
    Yc = m.Yc
    ym_c = Yc*(m.ξ*m.Pm_c*(1+m.τ_m_c)/m.ω_m_c)^(-m.σ)
    ym_k = m.Yk*(m.Pk^m.σ)*(m.ξ*m.Pm_k*(1+m.τ_m_k)/m.ω_m_k)^(-m.σ)

    if s.tariffsincome_PE==0
        m.tariffsincome = m.τ_m_c*m.ξ*m.Pm_c*ym_c + m.τ_m_k*m.ξ*m.Pm_k*ym_k
    end

end
#Display
 if s.display==1
     show("------------------------------------------------")
     show("Guesses: Φ_h=$(m.Φ_h) , w= $(m.w) , r= $(m.r) , ξ= $(m.ξ) , Pk= $(m.Pk)")
 end


 #Fixed costs
 if s.fcost_fgoods==0 % If in units of labor
     m.F     = m.w*m.F_base;
  else
     m.F     = m.F_base;
 end

end
