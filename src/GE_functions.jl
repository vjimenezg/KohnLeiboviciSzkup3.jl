    #########################################
    # Static problem for firm/entrepreneur
    #########################################

    # Specification:
    # One production function with materials
    # No working capital constraint
    # Derivations available at Model_Dec2018.tex


    function KLS3_staticproblem(m,s,r)

    ## Useful objects

        #Assets
        r=merge((a_grid_mat = r.a_grid'*ones(1,s.z_grid_size),

        #Productivity
        z_grid_mat = r.z),r)

        #ϕ_h
        ϕ_h = m.ϕ_h #Guessed value

        #Objects common for X and NX
        cap_gain = (1-m.Pk/m.Pk_lag)
        const_σ = ((m.σ-1)/m.σ)^m.σ
        m=merge((rtilde_u = ((m.r+m.δ)+(1-m.δ)*cap_gain)*m.Pk_lag,),m)

        if m.θ<1+m.r
            r=merge((k_const = (1/m.Pk_lag)*((1+m.r)/(1+m.r-m.θ))*r.a_grid_mat,),r)
        else
            r=merge((k_const = ones(size(r.a_grid_mat)),),r)
        end

    ## Exporters
        ϕ_x = ϕ_h + (m.τ^(1-m.σ))*m.Yf*((m.ξ/(1+m.τ_x))^m.σ)

    #Unconstrained

        # marginal cost
        r=merge((μ_u = (1./r.z_grid_mat) .* ((m.Pk/m.α_m)^m.α_m) .* ((m.w/((1-m.α)*(1-m.α_m))).^((1-m.α)*(1-m.α_m))) .* ((m.rtilde_u./(m.α*(1-m.α_m))).^(m.α*(1-m.α_m))),

        k_x_u = (m.α*(1-m.α_m)./m.rtilde_u).*const_σ*ϕ_x*(r.μ_u.^(1-m.σ)),
        n_x_u = ((1-m.α)*(1-m.α_m)/m.w)*const_σ*ϕ_x*(r.μ_u.^(1-m.σ)),
        m_x_u = (m.α_m/m.Pk)*const_σ*ϕ_x*(r.μ_u.^(1-m.σ)),

        yd_x_u = const_σ*ϕ_h*(r.μ_u.^(-m.σ)),
        yf_x_u = const_σ*(m.ξ^m.σ)*m.Yf*((1+m.τ_x)^(-m.σ))*(m.τ^(-m.sigma))*(r.μ_u.^(-m.σ)),

        pd_x_u = m.σ/(m.σ-1)*r.μ_u,
        pf_x_u = m.σ/(m.σ-1)*m.τ/m.ξ*r.μ_u,

        #Solution
        k_x = r.k_x_u,
        n_x = r.n_x_u,
        m_x = r.m_x_u,
        yd_x = r.yd_x_u,
        yf_x = r.yf_x_u,
        pd_x = r.pd_x_u,
        pf_x = r.pf_x_u,
        const_x = zeros(size(r.k_x_u))),r)

    #Constrained
        if m.θ<1+m.r
            r=merge((rtilde_x_c = (((1./r.k_const) .* m.α * (1-m.α_m) * const_σ * ϕ_x) .* ((r.μ_u .* m.rtilde_u.^(-m.α*(1-m.α_m))).^(1-m.σ))).^(1/(1-m.α*(1-m.α_m)*(1-m.σ))),

            mu_x_c = (1./r.z_grid_mat) .* ((m.Pk/m.α_m)^m.α_m) .* ((m.w/((1-m.α)*(1-m.α_m))).^((1-m.α)*(1-m.α_m))).* ((r.rtilde_x_c./(m.α*(1-m.α_m))).^(m.α*(1-m.α_m))),

            n_x_c = ((1-m.α)*(1-m.α_m)/m.w)*const_σ*ϕ_x*(r.μ_x_c.^(1-m.σ)),
            m_x_c = (m.α_m/m.Pk)*const_σ*ϕ_x*(r.μ_x_c.^(1-m.σ)),

            yd_x_c = const_σ*ϕ_h*(μ_x_c.^(-m.σ)),
            yf_x_c = const_σ*(m.ξ^m.σ)*m.Yf*((1+m.τ_x)^(-m.σ))*(m.τ^(-m.σ))*(r.μ_x_c.^(-m.σ)),

            pd_x_c = m.σ/(m.σ-1).*r.μ_x_c,
            pf_x_c = m.σ/(m.σ-1)*m.τ/m.ξ.*r.μ_x_c),r)

            #Solution
                const_x .= r.k_const .< r.k_x_u
                k_x[r.k_const .< r.k_x_u] .= r.k_const
                n_x[r.k_const .< r.k_x_u] .= r.n_x_c
                m_x[r.k_const .< r.k_x_u] .= r.m_x_c
                yd_x[r.k_const .< r.k_x_u] .= r.yd_x_c
                yf_x[r.k_const .< r.k_x_u] .= r.yf_x_c
                pd_x[r.k_const .< r.k_x_u] .= r.pd_x_c
                pf_x[r.k_const .< r.k_x_u] .= pf_x_c

        end

        r=merge((k_x = k_x,
        n_x = n_x,
        m_x = m_x,
        yd_x = yd_x,
        yf_x = yf_x,
        pd_x = pd_x,
        pf_x = pf_x,
        const_x = const_x),r)

    #Profits
        r=merge((π_x = r.pd_x.*r.yd_x + m.ξ*r.pf_x.*r.yf_x - ((m.r+m.δ)+(1-m.δ)*cap_gain)*m.Pk_lag*r.k_x - m.w*r.n_x - m.w*m.F_base - m.Pk*r.m_x,),r)

    ## Non-exporters

        ϕ_nx = ϕ_h

    #Unconstrained
        r=merge((μ_u = (1./r.z_grid_mat) * ((m.Pk/m.α_m)^m.α_m) * ((m.w/((1-m.α)*(1-m.α_m))).^((1-m.α)*(1-m.α_m))) * ((m.rtilde_u/(m.α*(1-m.α_m)))^(m.α*(1-m.α_m))),

        k_nx_u = (m.α*(1-m.α_m)/m.rtilde_u)*const_σ*ϕ_nx*(r.μ_u.^(1-m.σ)),
        n_nx_u = ((1-m.α)*(1-m.α_m)/m.w)*const_σ*ϕ_nx*(r.μ_u.^(1-m.σ)),
        m_nx_u = (m.α_m/m.Pk)*const_σ*ϕ_nx*(μ_u.^(1-m.σ)),

        yd_nx_u = const_σ*ϕ_h*(μ_u.^(-m.σ)),
        yf_nx_u = zeros(size(r.yd_nx_u)),

        pd_nx_u = m.σ/(m.σ-1)*r.μ_u,
        pf_nx_u = zeros(size(r.pd_nx_u))),r)

        #Solution
        k_nx = r.k_nx_u,
        n_nx = r.n_nx_u,
        m_nx = r.m_nx_u,
        yd_nx = r.yd_nx_u,
        yf_nx = r.yf_nx_u,
        pd_nx = r.pd_nx_u,
        pf_nx = r.pf_nx_u,
        const_nx = zeros(size(r.k_nx_u))

    #Constrained
        if m.θ<1+m.r
            r=merge((rtilde_nx_c = ( ((1./r.k_const) * m.α * (1-m.α_m) * const_σ * ϕ_nx) .* ((r.μ_u * m.rtilde_u.^(-m.α*(1-m.α_m))).^(1-m.σ)) ).^(1/(1-m.α*(1-m.α_m)*(1-m.σ))),

            μ_nx_c = (1./r.z_grid_mat) .* ((m.Pk/m.α_m)^m.α_m) .* ((m.w/((1-m.α)*(1-m.α_m))).^((1-m.α)*(1-m.α_m))) .* ((r.rtilde_nx_c./(m.α*(1-m.α_m))).^(m.α*(1-m.α_m))),

            n_nx_c = ((1-m.α)*(1-m.α_m)/m.w)*const_σ*ϕ_nx*(r.μ_nx_c.^(1-m.σ)),
            m_nx_c = (m.α_m/m.Pk)*const_σ*ϕ_nx*(r.μ_nx_c.^(1-m.σ)),

            yd_nx_c = const_σ*ϕ_h*(r.μ_nx_c.^(-m.σ)),
            yf_nx_c = zeros(size(r.yd_nx_c)),

            pd_nx_c = m.σ/(m.σ-1)*r.μ_nx_c,
            pf_nx_c = zeros(size(r.pd_nx_c))),r)

            #Solution
            const_nx .= r.k_const .< r.k_nx_u
            k_nx[r.k_const .< r.k_nx_u] .= r.k_const
            n_nx[r.k_const .< r.k_nx_u] .= r.n_nx_c
            m_nx[r.k_const .< r.k_nx_u] .= r.m_nx_c
            yd_nx[r.k_const .< r.k_nx_u] .= r.yd_nx_c
            yf_nx[r.k_const .< r.k_nx_u] .= r.yf_nx_c
            pd_nx[r.k_const .< r.k_nx_u] .= r.pd_nx_c
            pf_nx[r.k_const .< r.k_nx_u] .= r.pf_nx_c
        end

        r=merge((k_nx = k_nx,
        n_nx = n_nx,
        m_nx = m_nx,
        yd_nx = yd_nx,
        yf_nx = yf_nx,
        pd_nx = pd_nx,
        pf_nx = pf_nx,
        const_nx = const_nx),r)

    #Profits
        r=merge((π_nx = r.pd_nx.*r.yd_nx - ((m.r+m.δ)+(1-m.δ)*cap_gain)*m.Pk_lag*r.k_nx - m.w*r.n_nx - m.Pk*r.m_nx,),r)

    ## Export decision
        e .= r.pi_x .>= r.pi_nx;
        r=merge((e=e,),r)
    #    r.k_const = []; r.n_x_c = []; r.m_x_c = []; r.yd_x_c = [];r.yf_x_c = []; r.pd_x_c = [];  r.pf_x_c = [];
    #    r.k_x_u = [];   r.n_x_u = []; r.m_x_u = []; r.yd_x_u = []; r.yf_x_u = []; r.pd_x_u = []; r.pf_x_u = [];
    #    r.n_nx_c = []; r.m_nx_c = []; r.yd_nx_c = []; r.yf_nx_c = []; r.pd_nx_c = []; r.pf_nx_c = [];
    #    r.k_nx_u = []; r.n_nx_u = []; r.m_nx_u = []; r.yd_nx_u = []; r.yf_nx_u = []; r.pd_nx_u = []; r.pf_nx_u = [];
return r

    end


function KLS3_GE_par(x,m,s,r)

#### Aggregate prices and quantities ####

#Guessed prices

m=merge((w = exp(x[1]),
Φ_h = exp(x[2]),
ξ = exp(x[3]),
Pk = exp(x[4]),

Pk_lag = m.Pk),m)


### Solution

# Tariff income (initialized to 0)
if s.tariffsincome == 1
    m.Yk = exp(x[5])
    m.Yc = exp(x[2])
    m.Φ_h = (m.ω_h_c^m.σ)*m.Yc + ((m.ω_h_k*m.Pk)^m.σ)*m.Yk

    #Yc =  (m.Φ_h - ( m.ω_m_k*m.Pk)^m.σ*m.Yk)/( m.ω_m_c^m.σ)
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
 if s.fcost_fgoods==0 # If in units of labor
     m.F     = m.w*m.F_base;
  else
     m.F     = m.F_base;
 end

 #Static problems for tradable and nontradable sectors
     r =merge((KLS3_staticproblem(m,s,r),),r)


# PENDING (PERSISTENT)
# persistent guessV
# if isempty(guessV)
#    guessV=r.pi_nx;
# end
guessV=r.π_nx

#Dynamic problem and simulation (No sunk costs)
r=merge((KLS3_dynamicproblem(m,s,r,guessV),),r)
guessV = r.v

# PENDING (PERSISTENT)
# persistent guessM
# if isempty(guessM)
#    guessM=[r.z_pi'; zeros(length(r.a_grid)-1,length(r.z_grid))]; #Initialize
end

guessM=[r.z_π'; zeros(length(r.a_grid)-1,length(r.z_grid))]

end
