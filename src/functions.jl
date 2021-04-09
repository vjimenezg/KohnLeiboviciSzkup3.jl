    #########################################
    # Measure, Static Problem, Static Problem Period 2 (Alt Timing), Dynamic Problem, Simulation, GE Par (NLSolve), GE Par (Standard)#
    #########################################

############################# Measure ###########################

    function KLS3_measure(s,r,guessM)

    ## Computational objects
    numZ = length(r.z_grid)
    numA = length(r.a_grid)

    #Productivity
    #π_1 = r.z_π
    P = r.z_P

    ##
    # Mnew = zeros(numA,numZ)
    # Mnew[1,:] = π_1'
    Mnew=guessM

    # Only valid with value function iteration
    ap_ind = r.ap_ind

    iter = 0
    diff_M = 1


    while diff_M>s.eps_sm

        Mold = Mnew
        Mnew = zeros(numA,numZ)

        #for i=1:numA #Old asset state

        for j=1:numZ #Old productivity state
            PP=P[j,:]

            #aa_v=ap_ind[:,j]
            for i=1:numA #Old asset state
                #aa=ap_ind[i,j]
                #aa=aa_v[i]
                Mnew[ap_ind[i,j],:] = Mnew[ap_ind[i,j],:] + Mold[i,j]*PP
            end
        end

        diff_M = sum(abs.(Mnew-Mold))         # new criterion

        iter = iter + 1

    end

    measure = Mnew


    return measure
end


############################# Static Problem ########################

# One production function with materials
# No working capital constraint
# Derivations available at Model_Dec2018.tex

function KLS3_staticproblem(m,s,r)

### Useful objects

        #Assets
        r=merge(r,(a_grid_mat = r.a_grid'*ones(1,s.z_grid_size),

        #Productivity
        z_grid_mat = r.z))

        #ϕ_h

        ϕ_h = m.ϕ_h #Guessed value

        #Objects common for X and NX
        cap_gain = (1-m.Pk/m.Pk_lag)
        const_σ = ((m.σ-1)/m.σ)^m.σ
        m=merge(m,(rtilde_u = max.(((m.r+m.δ)+(1-m.δ)*cap_gain)*m.Pk_lag,0.00000001),))
        #display("rtilde_u=$(m.rtilde_u)")
        if m.θ<1+m.r
            r=merge(r,(k_const = (1/m.Pk_lag)*((1+m.r)/(1+m.r-m.θ))*r.a_grid_mat,))
        else
            r=merge(r,(k_const = ones(size(r.a_grid_mat)),))
        end

    ## Exporters
        ϕ_x = ϕ_h + (m.τ^(1-m.σ))*m.Yf*((m.ξ/(1+m.τ_x))^m.σ)

    #Unconstrained

        # marginal cost

        #r=merge(r,(μ_u = (1 ./r.z_grid_mat) .* (max.(m.Pk/m.α_m,0).^m.α_m) .* (max.(m.w/((1-m.α)*(1-m.α_m)),0).^((1-m.α)*(1-m.α_m))) .* (max.(m.rtilde_u./(m.α*(1-m.α_m)),0).^(m.α*(1-m.α_m))),))

        #no fix (max)

        r=merge(r,(μ_u = (1 ./r.z_grid_mat) .* ((m.Pk/m.α_m).^m.α_m) .* ((m.w/((1-m.α)*(1-m.α_m))).^((1-m.α)*(1-m.α_m))) .* ((m.rtilde_u./(m.α*(1-m.α_m))).^(m.α*(1-m.α_m))),))

        #μ_u_b = (1 ./r.z) .* ((Pk_b/m.α_m).^m.α_m) .* ((w_b/((1-m.α)*(1-m.α_m))).^((1-m.α)*(1-m.α_m))) .* ((m.rtilde_u./(m.α*(1-m.α_m))).^(m.α*(1-m.α_m)))

        r=merge(r,(k_x_u = (m.α*(1-m.α_m)./m.rtilde_u).*const_σ*ϕ_x*(r.μ_u.^(1-m.σ)),
        n_x_u = ((1-m.α)*(1-m.α_m)/m.w)*const_σ*ϕ_x*(r.μ_u.^(1-m.σ)),
        m_x_u = (m.α_m/m.Pk)*const_σ*ϕ_x*(r.μ_u.^(1-m.σ)),

        yd_x_u = const_σ*ϕ_h*(r.μ_u.^(-m.σ)),
        yf_x_u = const_σ*(m.ξ^m.σ)*m.Yf*((1+m.τ_x)^(-m.σ))*(m.τ^(-m.σ))*(r.μ_u.^(-m.σ)),

        pd_x_u = m.σ/(m.σ-1)*r.μ_u,
        pf_x_u = m.σ/(m.σ-1)*m.τ/m.ξ*r.μ_u))

        #Solution
        r=merge(r,(k_x = r.k_x_u,
        n_x = r.n_x_u,
        m_x = r.m_x_u,
        yd_x = r.yd_x_u,
        yf_x = r.yf_x_u,
        pd_x = r.pd_x_u,
        pf_x = r.pf_x_u,
        const_x = BitArray(undef,size(r.k_x_u))))

    #Constrained
        if m.θ<1+m.r
            #r=merge(r,(rtilde_x_c = max.(((1 ./r.k_const) .* m.α * (1-m.α_m) * const_σ * ϕ_x) .* ((r.μ_u .* max.(m.rtilde_u,0).^(-m.α*(1-m.α_m))).^(1-m.σ)),0).^(1/(1-m.α*(1-m.α_m)*(1-m.σ))),))

            #no fix(max)
            r=merge(r,(rtilde_x_c = (((1 ./r.k_const) .* m.α * (1-m.α_m) * const_σ * ϕ_x) .* ((r.μ_u .* (m.rtilde_u).^(-m.α*(1-m.α_m))).^(1-m.σ))).^(1/(1-m.α*(1-m.α_m)*(1-m.σ))),))

            r=merge(r,(μ_x_c = (1 ./r.z_grid_mat) .* ((m.Pk/m.α_m)^m.α_m) .* ((m.w/((1-m.α)*(1-m.α_m))).^((1-m.α)*(1-m.α_m))).* ((r.rtilde_x_c./(m.α*(1-m.α_m))).^(m.α*(1-m.α_m))),))

            r=merge(r,(n_x_c = ((1-m.α)*(1-m.α_m)/m.w)*const_σ*ϕ_x*(r.μ_x_c.^(1-m.σ)),
            m_x_c = (m.α_m/m.Pk)*const_σ*ϕ_x*(r.μ_x_c.^(1-m.σ)),

            yd_x_c = const_σ*ϕ_h*(r.μ_x_c.^(-m.σ)),
            yf_x_c = const_σ*(m.ξ^m.σ)*m.Yf*((1+m.τ_x)^(-m.σ))*(m.τ^(-m.σ))*(r.μ_x_c.^(-m.σ)),

            pd_x_c = m.σ/(m.σ-1).*r.μ_x_c,
            pf_x_c = m.σ/(m.σ-1)*m.τ/m.ξ.*r.μ_x_c))

            #Solution
                r.const_x .= r.k_const .< real(r.k_x_u)
                r.k_x[r.const_x] .= r.k_const[r.const_x]
                r.n_x[r.const_x] .= r.n_x_c[r.const_x]
                r.m_x[r.const_x] .= r.m_x_c[r.const_x]
                r.yd_x[r.const_x] .= r.yd_x_c[r.const_x]
                r.yf_x[r.const_x] .= r.yf_x_c[r.const_x]
                r.pd_x[r.const_x] .= r.pd_x_c[r.const_x]
                r.pf_x[r.const_x] .= r.pf_x_c[r.const_x]

        end


    #Profits
        r=merge(r,(π_x = r.pd_x.*r.yd_x .+ m.ξ*r.pf_x.*r.yf_x .- ((m.r+m.δ)+(1-m.δ)*cap_gain)*m.Pk_lag*r.k_x .- m.w*r.n_x .- m.w*m.F_base .- m.Pk*r.m_x,))

    ## Non-exporters

        ϕ_nx = ϕ_h

    #Unconstrained
        #r=merge(r,(μ_u = (1 ./r.z_grid_mat) * (max.(m.Pk/m.α_m,0)^m.α_m) * (max.(m.w/((1-m.α)*(1-m.α_m)),0).^((1-m.α)*(1-m.α_m))) * (max.(m.rtilde_u/(m.α*(1-m.α_m)),0).^(m.α*(1-m.α_m))),))

        #no fix(max)
        r=merge(r,(μ_u = (1 ./r.z_grid_mat) * ((m.Pk/m.α_m)^m.α_m) * ((m.w/((1-m.α)*(1-m.α_m))).^((1-m.α)*(1-m.α_m))) * ((m.rtilde_u/(m.α*(1-m.α_m))).^(m.α*(1-m.α_m))),))

        r=merge(r,(k_nx_u = (m.α*(1-m.α_m)/m.rtilde_u)*const_σ*ϕ_nx*(r.μ_u.^(1-m.σ)),
        n_nx_u = ((1-m.α)*(1-m.α_m)/m.w)*const_σ*ϕ_nx*(r.μ_u.^(1-m.σ)),
        m_nx_u = (m.α_m/m.Pk)*const_σ*ϕ_nx*(r.μ_u.^(1-m.σ)),

        yd_nx_u = const_σ*ϕ_h*(r.μ_u.^(-m.σ))))
        r=merge(r,(yf_nx_u = zeros(size(r.yd_nx_u)),

        pd_nx_u = m.σ/(m.σ-1)*r.μ_u))
        r=merge(r,(pf_nx_u = zeros(size(r.pd_nx_u)),))

        #Solution
        r=merge(r,(k_nx = r.k_nx_u,
        n_nx = r.n_nx_u,
        m_nx = r.m_nx_u,
        yd_nx = r.yd_nx_u,
        yf_nx = r.yf_nx_u,
        pd_nx = r.pd_nx_u,
        pf_nx = r.pf_nx_u,
        const_nx = BitArray(undef,size(r.k_nx_u))))

    #Constrained
        if m.θ<1+m.r
            #r=merge(r,(rtilde_nx_c = max.(((1 ./r.k_const) * m.α * (1-m.α_m) * const_σ * ϕ_nx) .* ((r.μ_u * max.(m.rtilde_u,0).^(-m.α*(1-m.α_m))).^(1-m.σ)),0).^(1/(1-m.α*(1-m.α_m)*(1-m.σ))),))

            #no fix (max)
            r=merge(r,(rtilde_nx_c = (((1 ./r.k_const) * m.α * (1-m.α_m) * const_σ * ϕ_nx) .* ((r.μ_u * (m.rtilde_u).^(-m.α*(1-m.α_m))).^(1-m.σ))).^(1/(1-m.α*(1-m.α_m)*(1-m.σ))),))

            r=merge(r,(μ_nx_c = (1 ./r.z_grid_mat) .* ((m.Pk/m.α_m)^m.α_m) .* ((m.w/((1-m.α)*(1-m.α_m))).^((1-m.α)*(1-m.α_m))) .* ((r.rtilde_nx_c./(m.α*(1-m.α_m))).^(m.α*(1-m.α_m))),))

            r=merge(r,(n_nx_c = ((1-m.α)*(1-m.α_m)/m.w)*const_σ*ϕ_nx*(r.μ_nx_c.^(1-m.σ)),
            m_nx_c = (m.α_m/m.Pk)*const_σ*ϕ_nx*(r.μ_nx_c.^(1-m.σ)),

            yd_nx_c = const_σ*ϕ_h*(r.μ_nx_c.^(-m.σ))))
            r=merge(r,(yf_nx_c = zeros(size(r.yd_nx_c)),

            pd_nx_c = m.σ/(m.σ-1)*r.μ_nx_c))
            r=merge(r,(pf_nx_c = zeros(size(r.pd_nx_c)),))

            #Solution
            r.const_nx .= r.k_const .< real(r.k_nx_u)
            r.k_nx[r.const_nx] .= r.k_const[r.const_nx]
            r.n_nx[r.const_nx] .= r.n_nx_c[r.const_nx]
            r.m_nx[r.const_nx] .= r.m_nx_c[r.const_nx]
            r.yd_nx[r.const_nx] .= r.yd_nx_c[r.const_nx]
            r.yf_nx[r.const_nx] .= r.yf_nx_c[r.const_nx]
            r.pd_nx[r.const_nx] .= r.pd_nx_c[r.const_nx]
            r.pf_nx[r.const_nx] .= r.pf_nx_c[r.const_nx]
        end



    #Profits
        r=merge(r,(π_nx = r.pd_nx.*r.yd_nx .- ((m.r+m.δ)+(1-m.δ)*cap_gain)*m.Pk_lag*r.k_nx .- m.w*r.n_nx .- m.Pk*r.m_nx,))

    ## Export decision
        r=merge(r,(e=BitArray(undef,size(r.k_nx_u)),))
        r.e .= real(r.π_x) .>= real(r.π_nx)

      #r=merge(r,(k_const = nothing, n_x_c = nothing, m_x_c = nothing, yd_x_c = nothing, yf_x_c = #nothing, pd_x_c = nothing,  pf_x_c = nothing))
      # r=merge(r,(k_x_u = nothing,   n_x_u = nothing, m_x_u = nothing, yd_x_u = nothing, yf_x_u = nothing, pd_x_u = nothing, pf_x_u = nothing))
      #r=merge(r,(n_nx_c = nothing, m_nx_c = nothing, yd_nx_c = nothing, yf_nx_c = nothing, pd_nx_c = nothing, pf_nx_c = nothing))
      # r=merge(r,(k_nx_u = nothing, n_nx_u = nothing, m_nx_u = nothing, yd_nx_u = nothing, yf_nx_u = nothing, pd_nx_u = nothing, pf_nx_u = nothing))

    return r

end
####################### Static Problem Period 2 #####################

function KLS3_staticproblem_period2_altTiming(m,s,r,rt)

    ### Useful objects

    #Assets
    r=merge(r,(a_grid_mat = r.a_grid'*ones(1,s.z_grid_size),

    #Productivity
    z_grid_mat = r.z))

    # capital gain
    cap_gain = (1-m.Pk/m.Pk_lag)

    # Capital
    r=merge(r,(k_x = rt[1].k,
    k_nx = rt[1].k))

    # Useful constants
    MC = (m.Pk/m.α_m )^m.α_m * (m.w /((1-m.α)*(1-m.α_m)))^((1-m.α)*(1-m.α_m))


### Exporters

    # useful constant
    const_μ_x = (1 ./r.z_grid_mat).^(1/(m.α.*(1-m.α_m))) .* 1 ./r.k_x .* ((m.σ-1)./(m.σ)).^m.σ .*(m.ϕ_h + (m.ξ/m.τ).^m.σ.*m.Yf.*m.τ.*(1+m.τ_x).^(-m.σ))

    # marginal cost
    r=merge(r,(μ_x = const_μ_x.^( (m.α.*(1-m.α_m)) ./ (1+(m.σ-1).*m.α.*(1-m.α_m)) ) .*MC.^(1 ./(1+(m.σ-1).*m.α.*(1-m.α_m))),))

    # variable inputs
    r=merge(r,(m_x = (m.α_m./m.Pk) .* (r.z_grid_mat.*r.μ_x).^(1 ./(m.α.*(1-m.α_m))) .* r.k_x .*MC.^(-1 ./(m.α.*(1-m.α_m))),
    n_x = ( ((1-m.α).*(1-m.α_m)) ./ m.w ) .* ( r.z_grid_mat.*r.μ_x ).^(1 ./(m.α.*(1-m.α_m))) .* r.k_x .*MC.^(-1 ./(m.α.*(1-m.α_m))),

    # output
    yd_x = ((m.σ-1)/m.σ).^m.σ .* m.ϕ_h .* r.μ_x.^(-m.σ),
    yf_x = ((m.σ-1)/m.σ).^m.σ .* m.ξ.^m.σ .* m.Yf .* (m.τ .*(1+m.τ_x) .*r.μ_x).^(-m.σ)))
    r=merge(r,(y_x = r.yd_x+m.τ*r.yf_x,

#Prices adjust (change in demand)
    pd_x = (r.yd_x./m.ϕ_h).^(-1/m.σ)))
    r.pd_x[r.yd_x.==0] .= 0
    r=merge(r,(pf_x = ((r.yf_x/m.Yf).^(-1/m.σ))/(1+m.τ_x),))
    r.pf_x[r.yf_x.==0] .= 0

#Profits (in units of the final good)
    r=merge(r,(π_x = r.pd_x.*r.yd_x .+ m.ξ*r.yf_x.*r.pf_x .- ((m.r+m.δ).+(1-m.δ)*cap_gain)*m.Pk_lag*r.k_x .- m.w*r.n_x .- m.Pk*r.m_x .- m.w*m.F_base,))


# Non-Exporters

    # useful constant
    const_μ_nx = (1 ./r.z_grid_mat).^(1/(m.α.*(1-m.α_m))) .* 1 ./r.k_nx .* ((m.σ-1)./(m.σ)).^m.σ .* m.ϕ_h

    # marginal cost
    r=merge(r,(μ_nx = const_μ_nx.^((m.α.*(1-m.α_m))./(1+(m.σ-1).*m.α.*(1-m.α_m)))*MC.^(1 ./(1 + (m.σ-1).*m.α.*(1-m.α_m))),))

    # variable inputs
    r=merge(r,(m_nx = (m.α_m/m.Pk) .* (r.z_grid_mat.*r.μ_nx).^(1 ./(m.α*(1-m.α_m))) .* r.k_nx .*MC.^(-1 ./(m.α.*(1-m.α_m))),
    n_nx = ( ((1-m.α).*(1-m.α_m)) ./ m.w ) .* ( r.z_grid_mat.*r.μ_nx ).^(1 ./(m.α.*(1-m.α_m))) .* r.k_nx .*MC.^(-1 ./(m.α.*(1-m.α_m))),

    # output
    yd_nx = ((m.σ-1)/m.σ).^m.σ .* m.ϕ_h .* r.μ_nx.^(-m.σ)))
    r=merge(r,(yf_nx = zeros(size(r.yd_nx)),))
    r=merge(r,(y_nx = (r.yd_nx+m.τ*r.yf_nx),

#Prices
    pd_nx = ((r.yd_nx/m.ϕ_h).^(-1/m.σ))))
	r.pd_nx[r.yd_nx.==0] .= 0
    r=merge(r,(pf_nx = zeros(s.a_grid_size,s.z_grid_size),

#Profits
    π_nx =  r.pd_nx.*r.yd_nx  .-  ((m.r+m.δ).+(1-m.δ)*cap_gain)*r.k_nx*m.Pk_lag .- m.w*r.n_nx .- m.Pk*r.m_nx,


# Export decision
    e=BitArray(undef,s.a_grid_size,s.z_grid_size)))
    r.e .= r.π_x .>= r.π_nx

    return r
end

########################### Dynamic Problem ########################

function KLS3_dynamicproblem(m,s,r,guessV)

     ### Dynamic problem: value functions and policy functions


     ###Initialize solution objects
        v_new = copy(guessV) #r.π_nx
        ap = zeros(Float64,size(v_new))
        ap_ind = zeros(Int64,size(v_new))
        c = zeros(Float64,size(v_new))


        v_old = copy(v_new)
        v_diff = 1
        z_P=r.z_P
        a_grid=r.a_grid
        exponentg=1-m.γ

        #Value function iteration algorithm

        profits = m.w .+ m.tariffsincome .+ r.e.*r.π_x + (1 .-r.e).*r.π_nx

        iter = 0
        while v_diff>s.eps

            iter +=1
            v_p=copy(v_old')

            for j = 1:s.z_grid_size

                v_pp = m.β*z_P[j,:]'*v_p
                for i = 1:s.a_grid_size

                    c1 =  profits[i,j] .+ a_grid[i].*(1 +m.r) .- a_grid
                    neg_c1_indexes = c1.<=0 #indices of a' for which consumption<0
                    u = (c1.^exponentg)./exponentg .+ neg_c1_indexes.*-1e50

                     v,index_ap = findmax(u .+ v_pp)

                      v_new[i,j] = v
                    ap_ind[i,j]=index_ap[2]

                end
            end
            for j =1:s.z_grid_size
                ap[:,j]=r.a_grid[ap_ind[:,j]]
            end

            # Accelerator

            # Consumption and utility
             c = profits .+ r.a_grid'.*ones(1,s.z_grid_size).*(1+m.r) .- ap  # The key is ap[:,:], the optimal asset choice

            u = ((c).^(1-m.γ))./(1-m.γ)

            for g=1:s.ac_iter
                for j = 1:s.z_grid_size
                    v_new[:,j] = u[:,j] .+ (m.β * z_P[j,:]'*v_new[ap_ind[:,j],:]')'

                end


            end

            # Convergence
             v_diff = maximum(abs.(v_new-v_old))
                v_old = copy(v_new)

            if mod(iter,100)==0
                show("Delta V: $v_diff")
            end

        end



    #Store output from value function iteration
        r=merge(r,(v = v_new,
        ap = ap))
        r=merge(r,(ap_ind=ap_ind,
        c=c,

    ## Store output to be used in simulation

        pd = (1 .-r.e).*r.pd_nx + r.e.*r.pd_x,
        yd = (1 .-r.e).*r.yd_nx + r.e.*r.yd_x,
        pf = r.e.*r.pf_x,
        yf = r.e.*r.yf_x,
        k = (1 .-r.e).*r.k_nx + r.e.*r.k_x,
        n = (1 .-r.e).*r.n_nx + r.e.*r.n_x,
        m = (1 .-r.e).*r.m_nx + r.e.*r.m_x,

        pd_nx = nothing, pd_x = nothing, yd_nx=nothing,  yd_x=nothing, pf_x=nothing, yf_x=nothing, k_nx =nothing, k_x = nothing, n_nx=nothing, n_x=nothing, m_nx=nothing, m_x=nothing,

        π_1 = (1 .-r.e).*r.π_nx + r.e.*r.π_x,
        a = r.a_grid_mat,

        #Fixed costs
        S_mat = r.e*m.F,
        F_mat = r.e*m.F))

    return r
end

 ############################## Simulation ##########################
function KLS3_simulate(m,s,r,guessM)


    # Compute stationary measure across states

     #guess=


#      if exist('sim.measure','var')
#
#          guessM=sim.measure
#
#      else
#
#         guessM=Mnew
#
#      end


    sim=(measure=KLS3_measure(s,r,guessM),

    w=m.w,
    ξ=m.ξ)

    # Exporters and non-exporters
    sim=merge(sim,(mass_x = sum(sim.measure.*r.e),))
    sim=merge(sim,(share_x=sim.mass_x, # share of exporters as a fraction of all entrepreneurs
    mass_d = sum(sim.measure.*(1 .-r.e)), #Mass of tradable producers who are non-exporters

    # Capital good

    I = m.δ*sum(sim.measure.*r.k),
    M = sum(sim.measure.*r.m)))

    yd_k = ((r.pd/(m.Pk*m.ω_h_k)).^(-m.σ)) * (sim.I + sim.M)
    ym_k = ((sim.ξ*m.Pm_k*(1+m.τ_m_k)/(m.Pk*m.ω_m_k)).^(-m.σ)) * (sim.I + sim.M)

    sim=merge(sim,(Yk = (sum(sim.measure.*m.ω_h_k.*(yd_k.^((m.σ-1)/m.σ)) ) + m.ω_m_k*(ym_k.^((m.σ-1)/m.σ)) )^(m.σ/(m.σ-1)),
    Pk = (sum(sim.measure.*(m.ω_h_k^m.σ).*(r.pd.^(1-m.σ))) + (m.ω_m_k^m.σ)*((sim.ξ*m.Pm_k*(1+m.τ_m_k)).^(1-m.σ)) )^(1/(1-m.σ)),

    # Consumption good

    C = sum(sim.measure.*r.c)))

    #yd_c = ((r.pd/m.ω_h_c).^(-m.σ)) * sim.C
    yd_c = max.(r.yd - yd_k,0.00001)
    ym_c = ((sim.ξ*m.Pm_c*(1+m.τ_m_c)/m.ω_m_c).^(-m.σ)) * sim.C

    sim=merge(sim,(Yc = (sum(sim.measure.*m.ω_h_c.*(yd_c.^((m.σ-1)/m.σ))) + m.ω_m_c*(ym_c^((m.σ-1)/m.σ)))^(m.σ/(m.σ-1)),)) #PENDING: FOR SOME REASON ym_c has negative values


    # ϕ_h

    sim=merge(sim,(ϕ_h = (m.ω_h_c^m.σ)*sim.Yc + ((m.ω_h_k*m.Pk)^m.σ)*sim.Yk,

    ### Compute price and quantity indexes

    #Domestic sales in units of final goods
     PdYd =sum(sim.measure[r.pd.>0].*r.pd[r.pd.>0].*r.yd[r.pd.>0]),

     #Exports
     PxYx = sim.ξ*sum(sim.measure[r.pf.>0].*r.pf[r.pf.>0].*r.yf[r.pf.>0])))
     sim=merge(sim,(PxYx_USD = sim.PxYx/sim.ξ,

     #Imports
    PmcYmc = sim.ξ*(1+m.τ_m_c)*m.Pm_c*ym_c,
    PmkYmk =  sim.ξ*(1+m.τ_m_k)*m.Pm_k*ym_k))
    sim=merge(sim,(PmYm = sim.PmcYmc + sim.PmkYmk,

    # Tariffs Income (T)
    tariffsincome=sim.ξ*m.τ_m_c*m.Pm_c*ym_c + sim.ξ*m.τ_m_k*m.Pm_k*ym_k,

    ### Compute aggregate variables needed to evaluate market clearing conditions

    #Labor and capital
    K = sum(sim.measure.*(r.k)))) #recall that already multiplied by oc_choice and sec_choice
    sim=merge(sim,(inv_agg = m.δ*sim.K,
    N = sum(sim.measure.*(r.n))))  # total labor demand for production

    # Fixed costs already corrected by sim.w if needed
    sim=merge(sim,(FC =  sum(sim.measure.*r.e.*r.F_mat),))

    # Total labor demand
    sim=merge(sim,(N = sim.N+((sim.FC)/sim.w)*(1-s.fcost_fgoods),


    ### Market clearing conditions

    #Labor
    n_supply = 1))
    sim=merge(sim,(n_demand = sim.N,))
    #sim=merge(sim,(mc_n = ((1+sim.n_demand)/(1+sim.n_supply))-1, #10*log.((1+sim.n_demand)/(1+sim.n_supply)),))
    sim=merge(sim,(mc_n = log.(sim.n_demand/sim.n_supply),

    #Assets
    #This market clearing condition is the result of reformulating the debt market clearing condition
    #In the original model, the sum of debt has to equal zero. Once the model is reformulated, this condition becomes the one below.
    a_supply = sum(sim.measure.*r.a_grid_mat),
    a_demand = sim.K))
    #sim=merge(sim,(mc_a = ((1+sim.a_demand)/(1+sim.a_supply)) - 1, #10*log.((1+sim.a_demand)/(1+sim.a_supply)),))
    sim=merge(sim,(mc_a = log.(sim.a_demand/sim.a_supply),

    #Final goods
    #mc_y = ((1+sim.y_demand)/(1+sim.y_supply))-1, #10*log.((1+sim.y_demand)/(1+sim.y_supply)),
    y_supply = sim.Yc,
    y_demand = sim.C)) #PENDING: Apparently has negative entries
    sim=merge(sim,(mc_y = log.(sim.y_demand/sim.y_supply),

    k_supply = sim.Yk,
    k_demand = sim.I + sim.M + sim.FC*s.fcost_fgoods))
    sim=merge(sim,(mc_k = log.(sim.k_demand/sim.k_supply),))



    #Beliefs
    #To solve the entrepreneur's problem, they need to have a belief about the aggregate price and quantity indexes in the market
    #In equilibrium, these beliefs need to be consistent with the actual aggregate prices and quantities
    #sim=merge(sim,(mc_y_belief = ((1+sim.ϕ_h)/(1+m.ϕ_h))-1,))
    #

    if s.tariffsincome==1
        sim=merge(sim,(mc_Yk_belief = log.(sim.Yk/m.Yk), #PENDING: apparently one of the Yks have negative entries
        mc_tariffs_belief = log.(sim.tariffsincome/m.tariffsincome),
        mc_y_belief = ((1+sim.Yc)/(1+m.Yc))-1))

    else
        sim=merge(sim,(mc_Yk_belief = 0,
        mc_tariffs_belief = 0,
        mc_y_belief = log.(sim.ϕ_h/m.ϕ_h)))
    end

    ## Main statistics

    sim=merge(sim,(Sales = sim.PdYd + sim.PxYx,)) #sim.PxYx is already denominated in units of the domestic final good
    sim=merge(sim,(GDP = sim.Sales - sim.M*m.Pk, # Value Added
    gross_output_ratio = (sim.ϕ_h)/(m.ξ*m.Yf),
    gross_output_ratio_2 = (m.Yc+m.Pk*m.Yk)/(m.ξ*m.Yf)))
    sim=merge(sim,(gross_output_ratio_3 =  sim.GDP/(sim.ξ*m.Yf),

    X_GDP = sim.PxYx/sim.GDP, # Exports, not value added exports
    D_GDP = sim.PdYd/sim.GDP,

    X_Sales = sim.PxYx/sim.Sales,
    D_Sales = sim.PdYd/sim.Sales,

    Absorption = sim.C+sim.I*sim.Pk,


    NX_GDP = (sim.PxYx-sim.PmYm)/sim.GDP)) #exports are pre-tax, imports include import tax
    sim=merge(sim,(NX_Absorption = (sim.PxYx-sim.PmYm)/sim.Absorption, #exports are pre-tax, imports include import tax
    NX_Sales = (sim.PxYx-sim.PmYm)/sim.Sales,

    X_D = sim.PxYx/sim.PdYd,

    IM_Sales = sim.PmYm/sim.Sales,
    IM_GDP = sim.PmYm/sim.GDP,
    IM_Absorption = sim.PmYm/sim.Absorption,

    Cimp_share = sim.PmcYmc/sim.PmYm,
    Kimp_share = sim.PmkYmk/sim.PmYm,

    Kimp_PkYk = sim.PmkYmk/(sim.Yk*sim.Pk),

    Kimp_Sales = sim.PmkYmk/sim.Sales,
    Kimp_GDP = sim.PmkYmk/sim.GDP,

    CRatio = sim.C/(sim.C+sim.I*sim.Pk),
    InvGDP=sim.inv_agg*sim.Pk/sim.GDP))

    # Every credit statistic
    r=merge(r,(d = (1+m.r)*(m.Pk*r.k .- r.a),))
    sim=merge(sim,(credit = sum(sim.measure.*max.(r.d,0)),)) # only entrepreneurs/firms (for workers r.d is negative)
    sim=merge(sim,(credit_gdp = sim.credit/sim.GDP,
    d_agg = sum(sim.measure.*r.d))) # for both firms and workers
    sim=merge(sim,(NFA_GDP = -sim.d_agg/sim.GDP,

    credit_sales = sim.credit/(sim.PdYd + sim.PxYx),

    k_wagebill = m.Pk*sim.K/(sim.w*sim.n_supply)))

    sales_x =  sim.ξ*r.pf.*r.yf #r.pf is denominated in foreign currency, so we adjust it
    sales_d = r.pd.*r.yd
    sales = sales_d+sales_x
    x_share = sales_x./sales

    sim=merge(sim,(x_share_av =  sum(filter(!isnan,sim.measure.*r.e.* x_share))./ sum(filter(!isnan,sim.measure.*r.e)),))
    #sim.x_share_av[isnan.(sim.x_share_av)].=1  #Don't know the role of this, but it throws error in Julia given sim.x_share_av is a scalar


    SalesExporters = sum(sim.measure.*sales.*r.e)
    sim=merge(sim,(x_share_wav =  Base.sum(filter(!isnan,sim.measure.*r.e.*(sales./SalesExporters).*x_share))/ sum(sim.measure.*r.e.*(sales./SalesExporters)),))
    #sim.x_share_wav[isnan.(x_share_wav)].=1 #Don't know the role of this, but it throws error in Julia given sim.x_share_av is a scalar

    ln_sales = log.(sales)
    ln_sales_d = log.(sales_d)
    ln_sales_ind=map(!,isnan.(ln_sales))
    ln_sales_d_ind=map(!,isnan.(ln_sales_d))

    sim=merge(sim,(ln_sales_mean = sum(sim.measure[ln_sales_ind].* ln_sales[ln_sales_ind]),))
    sim=merge(sim,(ln_sales_sd = sqrt.(sum(sim.measure[ln_sales_ind] .* (ln_sales[ln_sales_ind] .- sim.ln_sales_mean).^2)),

    ln_sales_d_mean = sum(sim.measure[ln_sales_d_ind] .* ln_sales_d[ln_sales_d_ind])))
    sim=merge(sim,(ln_sales_d_sd = sqrt.(sum(sim.measure[ln_sales_d_ind] .* (ln_sales_d[ln_sales_d_ind] .- sim.ln_sales_d_mean).^2)),))

    sales_mean = sum(sim.measure.*sales)
    sim=merge(sim,(sales_sd = sqrt.(sum(sim.measure .* (sales .- sales_mean).^2 )),))

    sales_d_mean = sum(sim.measure .*sales_d)
    sim=merge(sim,(sales_d_sd = sqrt.( sum(sim.measure .* (sales_d .- sales_d_mean).^2 )),))

    sim=merge(sim,(sales_avg_nx = sum(sim.measure.*(1 .-r.e).*sales) / sum(sim.measure.*(1 .-r.e)),
    labor_avg_nx = sum(sim.measure.*(1 .-r.e).*(r.n )) / sum(sim.measure.*(1 .-r.e))))
    sim=merge(sim,(sales_d_avg_nx = sim.sales_avg_nx,

    sales_avg_x = sum(sim.measure.*r.e.*sales)/ sum(sim.measure.*r.e),
    labor_avg_x = ((sim.FC/sim.w)*(1-s.fcost_fgoods) + sum(sim.measure.*r.e.*(r.n))) / sum(sim.measure.*r.e),

    # These include fixed costs
    labor_tot = sim.N,  # total labor demand for production
    labor_tot_x = (sim.FC/sim.w)*(1-s.fcost_fgoods)  + sum(sim.measure.*r.e.*(r.n)),
    labor_tot_nx =  sum(sim.measure.*(1 .-r.e).*(r.n)),
    labor_tot_w = sum(sim.measure)))

    sim=merge(sim,(labor_x_share = sim.labor_tot_x/sim.labor_tot,
    labor_nx_share = sim.labor_tot_nx/sim.labor_tot,

    sales_d_avg_x = sum(sim.measure.*r.e.*sales_d)/ sum(sim.measure.*r.e),

    xpremium_sales = sim.sales_avg_x/sim.sales_avg_nx,
    xpremium_labor = sim.labor_avg_x/sim.labor_avg_nx))

    sim=merge(sim,(xpremium_sales_d = sim.sales_d_avg_x/sim.sales_d_avg_nx,

    ln_sales_sd_mean = sim.ln_sales_sd /  sim.ln_sales_mean,
    ln_sales_d_sd_mean = sim.ln_sales_d_sd /  sim.ln_sales_d_mean,
    sales_sd_mean = sim.sales_sd /  sales_mean,
    sales_d_sd_mean = sim.sales_d_sd /  sales_d_mean))

    labor_mean= sum(sim.measure.*(r.n+ (r.e.*r.F_mat./sim.w )*(1-s.fcost_fgoods)))
    labor_sd=sqrt.(sum(sim.measure.*( r.n+ (r.e.*r.F_mat./sim.w)*(1-s.fcost_fgoods) .- labor_mean ).^2))
    sim=merge(sim,(labor_sd_mean = labor_sd /  labor_mean,


   # Average productivity
    avg_productivity1 = ((sum(sim.measure.*(r.z_grid_mat).^(m.σ-1)))/sum(sim.measure)).^(1/(m.σ - 1)),
    avg_productivity2 =  sum(sim.measure.*sales.* r.z_grid_mat) / sum(sim.measure .* sales),
    avg_productivity3 =  sum(sim.measure.*sales_d.* r.z_grid_mat) / sum(sim.measure .* sales_d),



 ## Solution-related statistics

# All
 a_min_share = sum(sim.measure[1,:]),
 a_max_share = sum(sim.measure[end,:])))

 return sim, r, s

end

############################# GE_par (2 versions, first for using NLsolve [in-place function], second standard)#############################

####### For NLsolve #######

function KLS3_GE_par!(F,x,m,s,r)

#### Aggregate prices and quantities ####

#Guessed prices

m = merge(m,(w = exp(x[1]),
ϕ_h = exp(x[2]),
ξ = exp(x[3]),
Pk = exp(x[4])))
m=merge(m,(Pk_lag = m.Pk,))


### Solution

# Tariff income (initialized to 0)
if s.tariffsincome == 1
    m=merge(m,(Yk = exp(x[5]),
    Yc = exp(x[2])))
    m=merge(m,(ϕ_h = (m.ω_h_c^m.σ)*m.Yc + ((m.ω_h_k*m.Pk)^m.σ)*m.Yk,))

    #Yc =  (m.ϕ_h - ( m.ω_m_k*m.Pk)^m.σ*m.Yk)/( m.ω_m_c^m.σ)
    Yc = m.Yc
    ym_c = Yc*(m.ξ*m.Pm_c*(1+m.τ_m_c)/m.ω_m_c)^(-m.σ)
    ym_k = m.Yk*(m.Pk^m.σ)*(m.ξ*m.Pm_k*(1+m.τ_m_k)/m.ω_m_k)^(-m.σ)

    if s.tariffsincome_PE==0
        m=merge(m,(tariffsincome = m.τ_m_c*m.ξ*m.Pm_c*ym_c + m.τ_m_k*m.ξ*m.Pm_k*ym_k,))
    end

end

#Display
 if s.display==1
     print("\nGuesses: ϕ_h=$(m.ϕ_h) , w= $(m.w) , r= $(m.r) , ξ= $(m.ξ) , Pk= $(m.Pk)")
 end


 #Fixed costs
 if s.fcost_fgoods==0 # If in units of labor
     m=merge(m,(F= m.w*m.F_base,))
  else
     m=merge(m,(F= m.F_base,))
 end

 #Static problems for tradable and nontradable sectors
     r_static=KLS3_staticproblem(m,s,r)
     r=merge(r,r_static)

# PENDING (PERSISTENT)
# persistent guessV
# if isempty(guessV)
#    guessV=r.π_nx
# end
guessV=copy(r.π_nx)

#Dynamic problem and simulation (No sunk costs)
r_dynamic=KLS3_dynamicproblem(m,s,r,guessV)
r=merge(r,r_dynamic)
guessV = r.v

# PENDING (PERSISTENT)
# persistent guessM
# if isempty(guessM)
#    guessM=vcat(r.z_π',zeros(s.a_grid_size-1,s.z_grid_size)) #Initialize
# end

guessM=vcat(r.z_π',zeros(s.a_grid_size-1,s.z_grid_size))


sim,r1,s1 = KLS3_simulate(m,s,r,guessM)
r=merge(r,r1)
s=merge(s,s1)

guessM=sim.measure

#Market clearing conditions
    if s.tariffsincome==0
        F[1]=sim.mc_n
        F[2]=sim.mc_y
        F[3]=sim.mc_y_belief
        F[4]=sim.mc_k
        sim=merge(sim,(mcc = F,))
    elseif s.tariffsincome==1
        F[1]=sim.mc_n
        F[2]=sim.mc_y
        F[3]=sim.mc_y_belief
        F[4]=sim.mc_k
        F[5]=sim.mc_Yk_belief
        #F[1]=sim.mc_n
        #F[2]=sim.mc_y
        #F[3]=sim.mc_y_belief
        #F[4]=sim.mc_k
        #F[5]=sim.mc_tariffs_belief
    end

#Display
if s.display==1

    print("GE: Y_MCC= $(sim.mc_y), N_MCC= $(sim.mc_n), A_MCC= $(sim.mc_a), Y_Belief= $(sim.mc_y_belief), K_MCC= $(sim.mc_k), Yk_MCC= $(sim.mc_Yk_belief)")
    print("\n")
    print("\nShare of exporters (all firms):  $(sim.share_x)")
    print("\nExporter domestic sales premium: $(sim.xpremium_sales_d)")
    print("\nCredit/GDP: $(sim.credit_gdp)")
    print("\nM/GDP: $(sim.IM_GDP)")
    print("\nX-M/GDP:  $(sim.NX_GDP)")
    print("\nPmcYmc/PmYm: $(sim.Cimp_share)")
#         print("\nM/Absorption:  $(sim.IM_Absorption )")
#         print("\nX-M/Absorption:  $(sim.NX_Absorption)")
#         print("\nM/Sales:  $(sim.IM_Sales)")
#         print("\nX-M/Sales:  $(sim.NX_Sales)")
#         print("\nPmkYmk/PkYk:  $(sim.Kimp_PkYk )")
#         print("\nPmkYmk/GDP:  $(sim.Kimp_GDP )")
#         print("\nPmkYmk/Sales:  $(sim.Kimp_Sales )")
#         print("\nPmcYmc/PmYm:  $(sim.Cimp_share )")
#         print("\nPmkYmk/PmYm:  $(sim.Kimp_share )")
#         print("\nX/GDP:  $(sim.X_GDP)")
#         print("\nX/Sales:  $(sim.X_Sales)")
#         print("\nAv. X intensity:  $(sim.x_share_av)")
#         print("\nC/(C+I):  $(sim.CRatio )")
#         print("\nInv_GDP:  $(sim.InvGDP )")
    print("\nShare of agents at highest a: $( sim.a_max_share)")

 end

    return F

end
####### Standard #######

function KLS3_GE_par(x,m,s,r)

#### Aggregate prices and quantities ####

#Guessed prices

m=merge(m,(w = exp(x[1]),
ϕ_h = exp(x[2]),
ξ = exp(x[3]),
Pk = exp(x[4])))

m=merge(m,(Pk_lag = m.Pk,))


### Solution

# Tariff income (initialized to 0)
if s.tariffsincome == 1
    m=merge(m,(Yk = exp(x[5]),
    Yc = exp(x[2])))
    m=merge(m,(ϕ_h = (m.ω_h_c^m.σ)*m.Yc + ((m.ω_h_k*m.Pk)^m.σ)*m.Yk,))

    #Yc =  (m.ϕ_h - ( m.ω_m_k*m.Pk)^m.σ*m.Yk)/( m.ω_m_c^m.σ)
    Yc = m.Yc
    ym_c = Yc*(m.ξ*m.Pm_c*(1+m.τ_m_c)/m.ω_m_c)^(-m.σ)
    ym_k = m.Yk*(m.Pk^m.σ)*(m.ξ*m.Pm_k*(1+m.τ_m_k)/m.ω_m_k)^(-m.σ)

    if s.tariffsincome_PE==0
        m=merge(m,(tariffsincome = m.τ_m_c*m.ξ*m.Pm_c*ym_c + m.τ_m_k*m.ξ*m.Pm_k*ym_k,))
    end

end

#Display
 if s.display==1
     print("\nGuesses: ϕ_h=$(m.ϕ_h) , w= $(m.w) , r= $(m.r) , ξ= $(m.ξ) , Pk= $(m.Pk)")
 end


 #Fixed costs
 if s.fcost_fgoods==0 # If in units of labor
     m=merge(m,(F= m.w*m.F_base,))
  else
     m=merge(m,(F = m.F_base,))
 end

 #Static problems for tradable and nontradable sectors
     r_static=KLS3_staticproblem(m,s,r)
     r=merge(r,r_static)

# PENDING (PERSISTENT)
# persistent guessV
# if isempty(guessV)
#    guessV=r.π_nx
# end
guessV=copy(r.π_nx)

#Dynamic problem and simulation (No sunk costs)
r_dynamic=KLS3_dynamicproblem(m,s,r,guessV)
r=merge(r,r_dynamic)
guessV = copy(r.v)

# PENDING (PERSISTENT)
# persistent guessM
# if isempty(guessM)
#    guessM=vcat(r.z_π',zeros(s.a_grid_size-1,s.z_grid_size)) #Initialize
# end

guessM=vcat(r.z_π',zeros(s.a_grid_size-1,s.z_grid_size))

sim,r1,s1 = KLS3_simulate(m,s,r,guessM)
r=merge(r,r1)
s=merge(s,s1)

guessM=sim.measure

#Market clearing conditions
    if s.tariffsincome==0
        mcc = [sim.mc_n sim.mc_y sim.mc_y_belief sim.mc_k]
        sim=merge(sim,(mcc = mcc,))
    elseif s.tariffsincome==1
        mcc = [sim.mc_n sim.mc_y sim.mc_y_belief sim.mc_k sim.mc_Yk_belief]
        # mcc = [sim.mc_n sim.mc_y sim.mc_y_belief sim.mc_k sim.mc_tariffs_belief]
        sim=merge(sim,(mcc = mcc,))
    end

#Display
if s.display==1

        print("\nGE: Y_MCC= $(sim.mc_y), N_MCC= $(sim.mc_n), A_MCC= $(sim.mc_a), Y_Belief= $(sim.mc_y_belief), K_MCC= $(sim.mc_k), Yk_MCC= $(sim.mc_Yk_belief)")
        print("\n")
        print("\nShare of exporters (all firms):  $(sim.share_x)")
        print("\nExporter domestic sales premium: $(sim.xpremium_sales_d)")
        print("\nCredit/GDP: $(sim.credit_gdp)")
        print("\nM/GDP: $(sim.IM_GDP)")
        print("\nX-M/GDP:  $(sim.NX_GDP)")
        print("\nPmcYmc/PmYm: $(sim.Cimp_share)")
#         print("\nM/Absorption:  $(sim.IM_Absorption )")
#         print("\nX-M/Absorption:  $(sim.NX_Absorption)")
#         print("\nM/Sales:  $(sim.IM_Sales)")
#         print("\nX-M/Sales:  $(sim.NX_Sales)")
#         print("\nPmkYmk/PkYk:  $(sim.Kimp_PkYk )")
#         print("\nPmkYmk/GDP:  $(sim.Kimp_GDP )")
#         print("\nPmkYmk/Sales:  $(sim.Kimp_Sales )")
#         print("\nPmcYmc/PmYm:  $(sim.Cimp_share )")
#         print("\nPmkYmk/PmYm:  $(sim.Kimp_share )")
#         print("\nX/GDP:  $(sim.X_GDP)")
#         print("\nX/Sales:  $(sim.X_Sales)")
#         print("\nAv. X intensity:  $(sim.x_share_av)")
#         print("\nC/(C+I):  $(sim.CRatio )")
#         print("\nInv_GDP:  $(sim.InvGDP )")
        print("\nShare of agents at highest a: $( sim.a_max_share)")

 end

 return mcc, m, r, s, sim

end
