############################# Functions used in the computation of the Transition (in addition to the ones in functions.jl) [measure_trans, dynamicproblem_trans_vec_t, static problem period 2 (alt timing), simulate_trans, transition_vec2 (NLSolve), transition_vec2 (standard)] ########################

########################## Measure Trans #############################
function KLS3_measure_trans(p,rt,N,M0)
@unpack z_grid,a_grid,z_P=p #allows to use the parameters in p without having to addend "p.". Furthermore, it makes clear which parameters are used in each function.

    numZ = length(z_grid)
    numA = length(a_grid)


    M = zeros(numA,numZ,N)
    M[:,:,1] = M0
    M[:,:,2] = M0 # in the second period news arrives after states are determined
    Mnew = copy(M0)

    P=z_P
    for n=2:N-1
        apt_ind = rt[n].ap_ind
        Mold = copy(Mnew)
        Mnew = zeros(numA,numZ)
        for j=1:numZ #Old productivity state
            PP=@view P[j:j,:] #@view (and @views) lets P[j:j,:] reference the variable P instead of allocating a new array, P[j:j,:]. This speeds up the code and is thus used several times.
            for i=1:numA #Old asset state
                @views Mnew[apt_ind[i,j]:apt_ind[i,j],:] += Mold[i,j]*PP
            end
        end

        M[:,:,n+1] = copy(Mnew)

    end
return M

end

####################### Static Problem Period 2 #####################

function KLS3_staticproblem_period2_altTiming(p_o_t,p)

@unpack a_grid,z_grid_size,z,α,α_m,τ,Yf,τ_x,σ,δ,F_base,a_grid_size,r = p #See comment at the beginning of KL3_measure_trans
    ### Useful objects

    #Assets
    p_o_t=merge(p_o_t,(a_grid_mat = a_grid'*ones(1,z_grid_size),

    #Productivity
    z_grid_mat = z))

    # capital gain
    cap_gain = (1-p_o_t.Pk/p_o_t.Pk_lag)

    # Capital
    p_o_t=merge(p_o_t,(k_x = rt[1].k,
    k_nx = rt[1].k))

    # Useful constants
    MC = (p_o_t.Pk/α_m )^α_m * (p_o_t.w /((1-α)*(1-α_m)))^((1-α)*(1-α_m))


### Exporters

    # useful constant
    const_μ_x = (1 ./p_o_t.z_grid_mat).^(1/(α.*(1-α_m))) .* 1 ./p_o_t.k_x .* ((σ-1)./(σ)).^σ .*(p_o_t.ϕ_h + (p_o_t.ξ/τ).^σ.*Yf.*τ.*(1+τ_x).^(-σ))

    # marginal cost
    p_o_t=merge(p_o_t,(μ_x = const_μ_x.^( (α.*(1-α_m)) ./ (1+(σ-1).*α.*(1-α_m)) ) .*MC.^(1 ./(1+(σ-1).*α.*(1-α_m))),))

    # variable inputs
    p_o_t=merge(p_o_t,(m_x = (α_m./p_o_t.Pk) .* (p_o_t.z_grid_mat.*p_o_t.μ_x).^(1 ./(α.*(1-α_m))) .* p_o_t.k_x .*MC.^(-1 ./(α.*(1-α_m))),
    n_x = ( ((1-α).*(1-α_m)) ./ p_o_t.w ) .* ( p_o_t.z_grid_mat.*p_o_t.μ_x ).^(1 ./(α.*(1-α_m))) .* p_o_t.k_x .*MC.^(-1 ./(α.*(1-α_m))),

    # output
    yd_x = ((σ-1)/σ).^σ .* p_o_t.ϕ_h .* p_o_t.μ_x.^(-σ),
    yf_x = ((σ-1)/σ).^σ .* p_o_t.ξ.^σ .* Yf .* (τ .*(1+τ_x) .*p_o_t.μ_x).^(-σ)))
    p_o_t=merge(p_o_t,(y_x = p_o_t.yd_x+τ*p_o_t.yf_x,

#Prices adjust (change in demand)
    pd_x = (p_o_t.yd_x./p_o_t.ϕ_h).^(-1/σ)))
    p_o_t.pd_x[p_o_t.yd_x.==0] .= 0
    p_o_t=merge(p_o_t,(pf_x = ((p_o_t.yf_x/Yf).^(-1/σ))/(1+τ_x),))
    p_o_t.pf_x[p_o_t.yf_x.==0] .= 0

#Profits (in units of the final good)
    p_o_t=merge(p_o_t,(π_x = p_o_t.pd_x.*p_o_t.yd_x .+ p_o_t.ξ*p_o_t.yf_x.*p_o_t.pf_x .- ((r+δ).+(1-δ)*cap_gain)*p_o_t.Pk_lag*p_o_t.k_x .- p_o_t.w*p_o_t.n_x .- p_o_t.Pk*p_o_t.m_x .- p_o_t.w*F_base,))


# Non-Exporters

    # useful constant
    const_μ_nx = (1 ./p_o_t.z_grid_mat).^(1/(α.*(1-α_m))) .* 1 ./p_o_t.k_nx .* ((σ-1)./(σ)).^σ .* p_o_t.ϕ_h

    # marginal cost
    p_o_t=merge(p_o_t,(μ_nx = const_μ_nx.^((α.*(1-α_m))./(1+(σ-1).*α.*(1-α_m)))*MC.^(1 ./(1 + (σ-1).*α.*(1-α_m))),))

    # variable inputs
    p_o_t=merge(p_o_t,(m_nx = (α_m/p_o_t.Pk) .* (p_o_t.z_grid_mat.*p_o_t.μ_nx).^(1 ./(α*(1-α_m))) .* p_o_t.k_nx .*MC.^(-1 ./(α.*(1-α_m))),
    n_nx = ( ((1-α).*(1-α_m)) ./ p_o_t.w ) .* ( p_o_t.z_grid_mat.*p_o_t.μ_nx ).^(1 ./(α.*(1-α_m))) .* p_o_t.k_nx .*MC.^(-1 ./(α.*(1-α_m))),

    # output
    yd_nx = ((σ-1)/σ).^σ .* p_o_t.ϕ_h .* p_o_t.μ_nx.^(-σ)))
    p_o_t=merge(p_o_t,(yf_nx = zeros(size(p_o_t.yd_nx)),))
    p_o_t=merge(p_o_t,(y_nx = (p_o_t.yd_nx+τ*p_o_t.yf_nx),

#Prices
    pd_nx = ((p_o_t.yd_nx/p_o_t.ϕ_h).^(-1/σ))))
	p_o_t.pd_nx[p_o_t.yd_nx.==0] .= 0
    p_o_t=merge(p_o_t,(pf_nx = zeros(a_grid_size,z_grid_size),

#Profits
    π_nx =  p_o_t.pd_nx.*p_o_t.yd_nx  .-  ((r+δ).+(1-δ)*cap_gain)*p_o_t.k_nx*p_o_t.Pk_lag .- p_o_t.w*p_o_t.n_nx .- p_o_t.Pk*p_o_t.m_nx,


# Export decision
    e=BitArray(undef,a_grid_size,z_grid_size))) #BitArray: Array of only 1's and 0's, uses less memory than standard Float64 arrays (hence, makes the code faster)
    p_o_t.e .= p_o_t.π_x .>= p_o_t.π_nx

    return p_o_t
end


############## Dynamic Problem Transition Vectorized t ###############

function KLS3_dynamicproblem_trans_vec_t(p_o_t,p,rt)

 @unpack N,a_grid_size,z_grid_size,a_grid,r,γ,z_P,a_grid_vec,β = p #See comment at the beginning of KL3_measure_trans

### Dynamic problem: value functions and policy functions

### Shocks

#Generate solution objects
    v_p=copy(rt[N].v)
    if s.transition_AD == 1 # See comment in line 202 of main.jl.
    v_new = zeros(eltype(p_o_t.w),size(v_p))
    else
    v_new = zeros(size(v_p))
    end

    ap_ind = zeros(Int64,size(v_new))
    ap=zeros(a_grid_size,z_grid_size)
    asset_income = a_grid'.*(1+r)
    exponentg=1-γ
### Dynamic problem: value functions and policy functions
    P=z_P
    v_pt=v_p'
    ones_vec=β*ones(a_grid_size,1)
    neg_c_indexes=BitArray(undef,a_grid_size,z_grid_size) #BitArray: Array of only 1's and 0's, uses less memory than standard Float64 arrays (hence, makes the code faster)

    #Value function iteration algorithm

    for t=N-1:-1:2
        profits = rt[t].profits
        for j = 1:z_grid_size

                c =  profits[:,j] + asset_income .- a_grid_vec
                neg_c_indexes = real(c).<=0 #indices of a' for which consumption<0
                u = (c.^exponentg)/exponentg + neg_c_indexes*-1e50
                MM=ones_vec*P[j,:]'
                v,index_ap = findmax(real(u + MM*v_pt),dims=2)
                v_new[:,j] = v
                for (i,index) in enumerate(index_ap)
                ap_ind[i,j] = index[2]
                end
                ap[:,j]=a_grid[ap_ind[:,j]]
        end

        v_pt=copy(v_new') #Store output from value function iteration
        rt[t]=merge(rt[t],(v = copy(v_new),
        ap = copy(ap),
        ap_ind= copy(ap_ind)))
        rt[t]=merge(rt[t],(c = profits + rt[t].a_grid_mat.*(1+r) .-rt[t].ap,

    ## Store output to be used in simulation

        pd = (1 .-rt[t].e).*rt[t].pd_nx .+ rt[t].e.*rt[t].pd_x,
        yd = (1 .-rt[t].e).*rt[t].yd_nx .+ rt[t].e.*rt[t].yd_x,
        pf = rt[t].e.*rt[t].pf_x,
        yf = rt[t].e.*rt[t].yf_x,
        k = (1 .-rt[t].e).*rt[t].k_nx .+ rt[t].e.*rt[t].k_x,
        n_x = rt[t].e.*(rt[t].n_x),
        n_nx = (1 .-rt[t].e).*rt[t].n_nx))
        rt[t]=merge(rt[t],(n = (1 .-rt[t].e).*rt[t].n_nx .+ rt[t].e.*rt[t].n_x,
        m = (1 .-rt[t].e).*rt[t].m_nx .+ rt[t].e.*rt[t].m_x))

        rt[t]=merge(rt[t],(pd_nx = nothing,
        pd_x = nothing,
        yd_nx=nothing,
        yd_x=nothing,
        pf_x=nothing,
        yf_x=nothing,
        k_nx =nothing,
        k_x = nothing,
        n_nx=nothing,
        n_x=nothing,
        m_nx=nothing,
        m_x=nothing))

        rt[t]=merge(rt[t],(π_1 = (1 .-rt[t].e).*rt[t].π_nx .+ rt[t].e.*rt[t].π_x,
        a = copy(rt[t].a_grid_mat)))

    end

return rt
end
############## Simulate Transition ###############


function KLS3_simulate_trans(p,p_o_t,rt,Guess)
@unpack N,ω_h_c,σ,ω_h_k,rv,Pfv,Yfv,θ_v,β_v,δ_v,τ_v,τ_x_v,τ_m_c_v,τ_m_k_v,δ,Pm_k,τ_m_k,ω_m_k,Pm_c,τ_m_c,ω_m_c,r = p #See comment at the beginning of KL3_measure_trans

    wguess = [p_o_t.wt[1] exp.(Guess[N-1:2*N-4])' p_o_t.wt[N]]
    ξguess = [p_o_t.ξt[1] exp.(Guess[2*N-3:3*N-6])' p_o_t.ξt[N]]
    Pkguess = [p_o_t.Pkt[1] exp.(Guess[3*N-5:4*N-8])' p_o_t.Pkt[N]]

    if s.tariffsincome == 1
        Ycguess = [p_o_t.Yct[1] exp.(Guess[1:N-2])' p_o_t.Yct[N]]
        Ykguess = [p_o_t.Ykt[1] exp.(Guess[4*N-7:5*N-10])' p_o_t.Ykt[N]]
        ϕhguess= (ω_h_c.^σ).*Ycguess + ((ω_h_k.*Pkguess).^σ).*Ykguess

    else
        ϕhguess = [p_o_t.ϕht[1] exp.(Guess[1:N-2])' p_o_t.ϕht[N]]

    end


    #Measure of firms
    M0 = rt[1].measure
    sim=(measure = KLS3_measure_trans(p,rt,N,M0),)

    if s.transition_AD == 1 # See comment in line 202 of main.jl.
    sim=merge(sim,(w=zeros(eltype(p_o_t.w),N,1),ξ=zeros(eltype(p_o_t.w),N,1),Pk=zeros(eltype(p_o_t.w),N,1),ϕ_h=zeros(eltype(p_o_t.w),N,1),share_x=zeros(eltype(p_o_t.w),N,1),share_d=zeros(eltype(p_o_t.w),N,1),K=zeros(eltype(p_o_t.w),N,1),inv_agg=zeros(eltype(p_o_t.w),N,1),I=zeros(eltype(p_o_t.w),N,1),M=zeros(eltype(p_o_t.w),N,1),Yk=zeros(eltype(p_o_t.w),N,1),C=zeros(eltype(p_o_t.w),N,1),Yc=zeros(eltype(p_o_t.w),N,1),PdYd=zeros(eltype(p_o_t.w),N,1),PxYx=zeros(eltype(p_o_t.w),N,1),PxYx_USD=zeros(eltype(p_o_t.w),N,1),PmYm=zeros(eltype(p_o_t.w),N,1),tariffsincome=zeros(eltype(p_o_t.w),N,1),FC=zeros(eltype(p_o_t.w),N,1),N=zeros(eltype(p_o_t.w),N,1),n_supply=zeros(eltype(p_o_t.w),N,1),n_demand=zeros(eltype(p_o_t.w),N,1),mc_n=zeros(eltype(p_o_t.w),N,1),a_supply=zeros(eltype(p_o_t.w),N,1),a_demand=zeros(eltype(p_o_t.w),N,1),mc_a=zeros(eltype(p_o_t.w),N,1),y_supply=zeros(eltype(p_o_t.w),N,1),y_demand=zeros(eltype(p_o_t.w),N,1),mc_y=zeros(eltype(p_o_t.w),N,1),k_supply=zeros(eltype(p_o_t.w),N,1),k_demand=zeros(eltype(p_o_t.w),N,1),mc_k=zeros(eltype(p_o_t.w),N,1),mc_Yk_belief=zeros(eltype(p_o_t.w),N,1),mc_tariffs_belief=zeros(eltype(p_o_t.w),N,1),mc_y_belief=zeros(eltype(p_o_t.w),N,1),Sales=zeros(eltype(p_o_t.w),N,1),GDP=zeros(eltype(p_o_t.w),N,1),X_GDP=zeros(eltype(p_o_t.w),N,1),D_GDP=zeros(eltype(p_o_t.w),N,1),X_Sales=zeros(eltype(p_o_t.w),N,1),D_Sales=zeros(eltype(p_o_t.w),N,1),NX_GDP=zeros(eltype(p_o_t.w),N,1),NX_Sales=zeros(eltype(p_o_t.w),N,1),X_D=zeros(eltype(p_o_t.w),N,1),credit=zeros(eltype(p_o_t.w),N,1),credit_gdp=zeros(eltype(p_o_t.w),N,1),d_agg=zeros(eltype(p_o_t.w),N,1),NFA_GDP=zeros(eltype(p_o_t.w),N,1),k_wagebill=zeros(eltype(p_o_t.w),N,1),x_share_av=zeros(eltype(p_o_t.w),N,1),ln_sales_sd=zeros(eltype(p_o_t.w),N,1),ln_sales_d_sd=zeros(eltype(p_o_t.w),N,1),sales_sd=zeros(eltype(p_o_t.w),N,1),sales_d_sd=zeros(eltype(p_o_t.w),N,1),sales_avg_nx=zeros(eltype(p_o_t.w),N,1),labor_avg_nx=zeros(eltype(p_o_t.w),N,1),sales_d_avg_nx=zeros(eltype(p_o_t.w),N,1),sales_avg_x=zeros(eltype(p_o_t.w),N,1),labor_avg_x=zeros(eltype(p_o_t.w),N,1),labor_tot=zeros(eltype(p_o_t.w),N,1),labor_tot_x=zeros(eltype(p_o_t.w),N,1),labor_tot_nx=zeros(eltype(p_o_t.w),N,1),labor_tot_w=zeros(eltype(p_o_t.w),N,1),labor_x_share=zeros(eltype(p_o_t.w),N,1),labor_nx_share=zeros(eltype(p_o_t.w),N,1),sales_d_avg_x=zeros(eltype(p_o_t.w),N,1),xpremium_sales=zeros(eltype(p_o_t.w),N,1),xpremium_labor=zeros(eltype(p_o_t.w),N,1),xpremium_sales_d=zeros(eltype(p_o_t.w),N,1),ln_sales_sd_mean=zeros(eltype(p_o_t.w),N,1),ln_sales_d_sd_mean=zeros(eltype(p_o_t.w),N,1),sales_sd_mean=zeros(eltype(p_o_t.w),N,1),sales_d_sd_mean=zeros(eltype(p_o_t.w),N,1),labor_sd_mean=zeros(eltype(p_o_t.w),N,1),avg_productivity1=zeros(eltype(p_o_t.w),N,1),avg_productivity2=zeros(eltype(p_o_t.w),N,1),avg_productivity3=zeros(eltype(p_o_t.w),N,1),X_Laspeyres=zeros(eltype(p_o_t.w),N,1),D_Laspeyres=zeros(eltype(p_o_t.w),N,1),Sales_Laspeyres=zeros(eltype(p_o_t.w),N,1),GDP_Laspeyres=zeros(eltype(p_o_t.w),N,1),Imports=zeros(eltype(p_o_t.w),N,1),Imports_C=zeros(eltype(p_o_t.w),N,1),Imports_K=zeros(eltype(p_o_t.w),N,1),Imports_C_Laspeyres=zeros(eltype(p_o_t.w),N,1),Imports_K_Laspeyres=zeros(eltype(p_o_t.w),N,1),Imports_Laspeyres=zeros(eltype(p_o_t.w),N,1),a_min_share=zeros(eltype(p_o_t.w),N,1),a_max_share=zeros(eltype(p_o_t.w),N,1)))
    else
    sim=merge(sim,(w=zeros(N,1),ξ=zeros(N,1),Pk=zeros(N,1),ϕ_h=zeros(N,1),share_x=zeros(N,1),share_d=zeros(N,1),K=zeros(N,1),inv_agg=zeros(N,1),I=zeros(N,1),M=zeros(N,1),Yk=zeros(N,1),C=zeros(N,1),Yc=zeros(N,1),PdYd=zeros(N,1),PxYx=zeros(N,1),PxYx_USD=zeros(N,1),PmYm=zeros(N,1),tariffsincome=zeros(N,1),FC=zeros(N,1),N=zeros(N,1),n_supply=zeros(N,1),n_demand=zeros(N,1),mc_n=zeros(N,1),a_supply=zeros(N,1),a_demand=zeros(N,1),mc_a=zeros(N,1),y_supply=zeros(N,1),y_demand=zeros(N,1),mc_y=zeros(N,1),k_supply=zeros(N,1),k_demand=zeros(N,1),mc_k=zeros(N,1),mc_Yk_belief=zeros(N,1),mc_tariffs_belief=zeros(N,1),mc_y_belief=zeros(N,1),Sales=zeros(N,1),GDP=zeros(N,1),X_GDP=zeros(N,1),D_GDP=zeros(N,1),X_Sales=zeros(N,1),D_Sales=zeros(N,1),NX_GDP=zeros(N,1),NX_Sales=zeros(N,1),X_D=zeros(N,1),credit=zeros(N,1),credit_gdp=zeros(N,1),d_agg=zeros(N,1),NFA_GDP=zeros(N,1),k_wagebill=zeros(N,1),x_share_av=zeros(N,1),ln_sales_sd=zeros(N,1),ln_sales_d_sd=zeros(N,1),sales_sd=zeros(N,1),sales_d_sd=zeros(N,1),sales_avg_nx=zeros(N,1),labor_avg_nx=zeros(N,1),sales_d_avg_nx=zeros(N,1),sales_avg_x=zeros(N,1),labor_avg_x=zeros(N,1),labor_tot=zeros(N,1),labor_tot_x=zeros(N,1),labor_tot_nx=zeros(N,1),labor_tot_w=zeros(N,1),labor_x_share=zeros(N,1),labor_nx_share=zeros(N,1),sales_d_avg_x=zeros(N,1),xpremium_sales=zeros(N,1),xpremium_labor=zeros(N,1),xpremium_sales_d=zeros(N,1),ln_sales_sd_mean=zeros(N,1),ln_sales_d_sd_mean=zeros(N,1),sales_sd_mean=zeros(N,1),sales_d_sd_mean=zeros(N,1),labor_sd_mean=zeros(N,1),avg_productivity1=zeros(N,1),avg_productivity2=zeros(N,1),avg_productivity3=zeros(N,1),X_Laspeyres=zeros(N,1),D_Laspeyres=zeros(N,1),Sales_Laspeyres=zeros(N,1),GDP_Laspeyres=zeros(N,1),Imports=zeros(N,1),Imports_C=zeros(N,1),Imports_K=zeros(N,1),Imports_C_Laspeyres=zeros(N,1),Imports_K_Laspeyres=zeros(N,1),Imports_Laspeyres=zeros(N,1),a_min_share=zeros(N,1),a_max_share=zeros(N,1)))
    end

    for t=2:N
    ### Shocks

        p_o_t=merge(p_o_t,(r = rv[t], # shock to interest rate
        r_lag = rv[t-1],
        Pf_lag = Pfv[t-1], # shock to foreign price
        Pf = Pfv[t],
        Yf = Yfv[t], # shock to foreign demand
        θ = θ_v[t], # shock to collateral constraint
        β = β_v[t],    # shock to discount factor
        δ = δ_v[t], # shock to depreciation rate
        τ = τ_v[t], # shock to iceberg costs
        τ_x = τ_x_v[t], # shocks to tariffs
        τ_m_c = τ_m_c_v[t], # shocks to consumption tariffs
        τ_m_k = τ_m_k_v[t], # shocks to capital tariffs


        ## GE prices
        ϕ_h = ϕhguess[t],
        w = wguess[t],
        ξ = ξguess[t],
        Pk = Pkguess[t],
        Pk_lag = Pkguess[t-1]))

        if s.tariffsincome==1
            p_o_t=merge(p_o_t,(Yk = Ykguess[t],
            Yc = Ycguess[t],
            tariffsincome=p_o_t.tariffsincomet[t]))
        end

        ## Useful objects

        n = rt[t].n
        mat = rt[t].m  # intermediates
        measure = sim.measure[:,:,t]  # measure
        pd = rt[t].pd
        yd = rt[t].yd
        pf = rt[t].pf
        yf = rt[t].yf
        e = rt[t].e
        k = rt[t].k
        if t<N
            kp = rt[t+1].k
        end
        F_mat = rt[t].F_mat
        c = rt[t].c

        # Saving prices

        sim.w[t,1]=p_o_t.w
        sim.ξ[t,1]=p_o_t.ξ
        sim.Pk[t,1] = p_o_t.Pk
        sim.ϕ_h[t,1] = p_o_t.ϕ_h

        # Exporters and non-exporters
        sim.share_x[t,1] = sum(measure.*e) # share of exporters
        sim.share_d[t,1] = sum(measure.*(1 .-e))  #share of non-exporters

        ## Capital good
        sim.K[t,1] = sum(real(measure.*k))
        if t>1 && t<N
            measurep =  sim.measure[:,:,t+1]
            sim.inv_agg[t,1] = real(sum(measurep.*kp) - (1 .-δ_v[t])*sum(measure.*k))
        elseif t==N #For last period (t=N)
            # We assume that we reached steady state when t=N
            sim.inv_agg[t,1] = δ*sim.K[t,1]
        end

        sim.I[t,1] = max(sim.inv_agg[t,1],0) #Demand for new investment goods
        sim.M[t,1] = sum(real(measure.*mat))

        yd_k = ((pd/(p_o_t.Pk*ω_h_k)).^(-σ)) * (sim.I[t,1] + sim.M[t,1])
        ym_k = ((p_o_t.ξ*Pm_k*(1+τ_m_k)/(p_o_t.Pk*ω_m_k)).^(-σ)) * (sim.I[t,1] + sim.M[t,1])

        sim.Yk[t,1] = (sum(real(measure.*ω_h_k.*(yd_k.^((σ-1)/σ)))) + ω_m_k*(ym_k.^((σ-1)/σ)))^(σ/(σ-1))
        sim.Pk[t,1] = (sum(real(measure.*(ω_h_k^σ).*(pd.^(1-σ)))) + (ω_m_k^σ)*((p_o_t.ξ*Pm_k*(1+τ_m_k)).^(1-σ)) )^(1/(1-σ))

        ## Consumption good

        sim.C[t,1] = sum(real(measure.*c))

        #yd_c = ((pd/ω_h_c).^(-σ)) * sim.C[t,1]
        yd_c = max.(real(yd - yd_k),0.00001)
        ym_c = (p_o_t.ξ*Pm_c*(1+τ_m_c)/(ω_m_c)).^(-σ) * sim.C[t,1]
        sim.Yc[t,1] = (sum(real(measure.*ω_h_c.*(yd_c.^((σ-1)/σ)))) + ω_m_c*(ym_c.^((σ-1)/σ)) )^(σ/(σ-1))


        ## ϕ_h

        sim.ϕ_h[t,1] = (ω_h_c^σ)*sim.Yc[t,1] + ((ω_h_k*p_o_t.Pk)^σ)*sim.Yk[t,1]


        ## Compute price and quantity indexes

        #Domestic sales
         sim.PdYd[t,1] =sum(real(measure[real(pd).>0].*(pd[real(pd).>0].*yd[real(pd).>0])))

        #Exports
         sim.PxYx[t,1] = p_o_t.ξ*sum(real(measure[real(pf).>0].*(pf[real(pf).>0].*yf[real(pf).>0])))

         sim.PxYx_USD[t,1] = sim.PxYx[t,1]/p_o_t.ξ

         #Imports
        sim.PmYm[t,1] = p_o_t.ξ*(1+τ_m_c)*Pm_c*ym_c + p_o_t.ξ*(1+τ_m_k)*Pm_k*ym_k

        # Tariffs Income (T)
        sim.tariffsincome[t,1]=p_o_t.ξ*τ_m_c*Pm_c*ym_c + p_o_t.ξ*τ_m_k*Pm_k*ym_k
        if s.tariffsincome_PE==1
            sim.tariffsincome[t,1] = p_o_t.tariffsincome
        end
        ## Compute aggregate variables needed to evaluate market clearing conditions

        #Labor and capital
        sim.FC[t,1] =  sum(measure.*e.*F_mat)
        sim.N[t,1] = sum(real(measure.*n))    +((sim.FC[t,1])/p_o_t.w)*(1 .-s.fcost_fgoods)

        ## Market clearing conditions

        #Labor
        sim.n_supply[t,1] = 1
        sim.n_demand[t,1] = sim.N[t,1]
        #sim.mc_n[t,1] = ((1+sim.n_demand[t,1])/(1+sim.n_supply[t,1]))-1
        sim.mc_n[t,1] = log.(sim.n_demand[t,1]/sim.n_supply[t,1])

        #Assets
        #This market clearing condition is the result of reformulating the debt market clearing condition
        #In the original model, the sum of debt has to equal zero. Once the model is reformulated, this condition becomes the one below.
        sim.a_supply[t,1] = sum(real(measure.*rt[t].a))
        sim.a_demand[t,1] = sim.K[t,1]
        #sim.mc_a[t,1] = (1+sim.a_demand[t,1])/(1+sim.a_supply[t,1])-1

        sim.mc_a[t,1] = log.(sim.a_demand[t,1]/sim.a_supply[t,1])

        #Final goods
        sim.y_supply[t,1] = sim.Yc[t,1]
        sim.y_demand[t,1] = sim.C[t,1]
        sim.mc_y[t,1] = log.(sim.y_demand[t,1]/sim.y_supply[t,1])

        # Capital goods
        sim.k_supply[t,1] = sim.Yk[t,1]
        sim.k_demand[t,1] = sim.I[t,1] + sim.M[t,1] + sim.FC[t,1]*s.fcost_fgoods
        sim.mc_k[t,1] = log.(sim.k_demand[t,1]/sim.k_supply[t,1])
#       sim.mc_k[t,1] =sim.k_demand[t,1]/sim.k_supply[t,1]-1

        #Beliefs
        #To solve the entrepreneur's problem, they need to have a belief about the aggregate price and quantity indexes in the market
        #In equilibrium, these beliefs need to be consistent with the actual aggregate prices and quantities
        #sim.mc_y_belief[t,1] = log.(sim.Yh[t,1]/m.Yh) sim.mc_y_belief = 10*log.((1+sim.Ycpi_t)/(1+m.Yh))  sim.mc_y_belief[t,1] = ((1+sim.Yh[t,1])/(1+m.Yh))-1


        if s.tariffsincome==1
            sim.mc_Yk_belief[t,1] = log.(sim.Yk[t,1]/p_o_t.Yk)
            sim.mc_tariffs_belief[t,1] = log.(sim.tariffsincome[t,1]/p_o_t.tariffsincome)
            sim.mc_y_belief[t,1] = ((1+sim.Yc[t,1])/(1+p_o_t.Yc))-1
        else
            sim.mc_Yk_belief[t,1] = 0
            sim.mc_tariffs_belief[t,1] = 0
            sim.mc_y_belief[t,1] = log.(sim.ϕ_h[t,1]/ϕhguess[t])
        end




        ## Main statistics

        sim.Sales[t,1] = sim.PdYd[t,1] + sim.PxYx[t,1] #sim.PxYx is already denominated in units of the domestic final good
        sim.GDP[t,1] = sim.Sales[t,1] -  sim.M[t,1]*p_o_t.Pk #

        sim.X_GDP[t,1] = sim.PxYx[t,1]/sim.GDP[t,1]
        sim.D_GDP[t,1]  = sim.PdYd[t,1] /sim.GDP[t,1]

        sim.X_Sales[t,1] = sim.PxYx[t,1]/sim.Sales[t,1]
        sim.D_Sales[t,1] = sim.PdYd[t,1]/sim.Sales[t,1]

         sim.NX_GDP[t,1] = (sim.PxYx[t,1]-sim.PmYm[t,1])/sim.GDP[t,1]
         sim.NX_Sales[t,1] = (sim.PxYx[t,1]-sim.PmYm[t,1])/sim.Sales[t,1]

         sim.X_D[t,1] = sim.PxYx[t,1]/sim.PdYd[t,1]

        # Every credit statistic only for tradables
        rt[t]=merge(rt[t],(d = (1+r).*(p_o_t.Pk_lag*k - rt[t].a),))

        sim.credit[t,1] = sum(real(measure.*max.(real(rt[t].d),0))) # only entrepreneurs/firms (for workers r.d is negative)
        sim.credit_gdp[t,1] = sim.credit[t,1]/sim.GDP[t,1]
        sim.d_agg[t,1] = sum(real(measure.*rt[t].d)) # for both firms and workers
        sim.NFA_GDP[t,1] = -sim.d_agg[t,1]/sim.GDP[t,1]

        sim.k_wagebill[t,1] = p_o_t.Pk_lag*sim.K[t,1]/(p_o_t.w*sim.n_supply[t,1])

        # Sales
        sales_x =  p_o_t.ξ*rt[t].pf.*rt[t].yf #r.pf is denominated in foreign currency, so we adjust it
        sales_d = rt[t].pd.*rt[t].yd
        sales = sales_d+sales_x
        x_share = sales_x./sales

        sim.x_share_av[t,1] =  sum(filter(!isnan,real(measure.*e.* x_share))) / sum(real(measure.*e))
        # sim.x_share_av[isnan.(sim.x_share_av[t,1])].=1 #Don't know the role of this, but it throws error in Julia given sim.x_share_av is a scalar

        if s.extra_results==1
            ln_sales = log.(sales)
            ln_sales_d =log.(sales_d)
            ln_sales_ind = isnan.(ln_sales).==0
            ln_sales_d_ind = isnan.(ln_sales_d).==0

            ln_sales_mean = sum(real(measure[ln_sales_ind].*ln_sales[ln_sales_ind]))
            sim.ln_sales_sd[t,1] = sqrt.(sum(real(measure[ln_sales_ind] .* (ln_sales[ln_sales_ind] .- ln_sales_mean).^2)))


            ln_sales_d_mean = sum(real(measure[ln_sales_d_ind] .* ln_sales_d[ln_sales_d_ind]))
            sim.ln_sales_d_sd[t,1] = sqrt.( sum(real(measure[ln_sales_d_ind] .* (ln_sales_d[ln_sales_d_ind] .- ln_sales_d_mean).^2)))

            sales_mean = sum(real(measure.*sales))
            sim.sales_sd[t,1] = sqrt.(sum(real(measure .* (sales .- sales_mean).^2)))

            sales_d_mean = sum(real(measure .* sales_d))
            sim.sales_d_sd[t,1] = sqrt.(sum(real(measure .* (sales_d .- sales_d_mean).^2)))

            sim.sales_avg_nx[t,1] = sum(real(measure.*(1 .-e).*sales)) / sum(real(measure.*(1 .-e)))
            sim.labor_avg_nx[t,1] = sum(real(measure.*(1 .-e).*(n))) / sum(real(measure.*(1 .-e)))
            sim.sales_d_avg_nx[t,1] = sim.sales_avg_nx[t,1]

            sim.sales_avg_x[t,1] = sum(real(measure.*e.*sales))/ sum(real(measure.*e))
            sim.labor_avg_x[t,1] = ((sim.FC[t,1]/p_o_t.w)*(1-s.fcost_fgoods) + sum(real(measure.*e.*(n)))) / sum(real(measure.*e))

        # These include fixed costs
            sim.labor_tot[t,1] = sim.N[t,1]  # total labor demand for production
            sim.labor_tot_x[t,1] = (sim.FC[t,1]/p_o_t.w)*(1-s.fcost_fgoods)  + sum(real(measure.*e.*(n)))
            sim.labor_tot_nx[t,1] =  sum(real(measure.*(1 .-e).*(n)))
            sim.labor_tot_w[t,1] = sum(real(measure))

            sim.labor_x_share[t,1] = sim.labor_tot_x[t,1]/sim.labor_tot[t,1]
            sim.labor_nx_share[t,1] = sim.labor_tot_nx[t,1]/sim.labor_tot[t,1]

            sim.sales_d_avg_x[t,1] = sum(real(measure.*e.*sales_d))/ sum(real(measure.*e))

            sim.xpremium_sales[t,1] = sim.sales_avg_x[t,1]/sim.sales_avg_nx[t,1]
            sim.xpremium_labor[t,1] = sim.labor_avg_x[t,1]/sim.labor_avg_nx[t,1]

            sim.xpremium_sales_d[t,1] = sim.sales_d_avg_x[t,1]/sim.sales_d_avg_nx[t,1]

            sim.ln_sales_sd_mean[t,1] = sim.ln_sales_sd[t,1] /  ln_sales_mean
            sim.ln_sales_d_sd_mean[t,1] = sim.ln_sales_d_sd[t,1] /  ln_sales_d_mean
            sim.sales_sd_mean[t,1] = sim.sales_sd[t,1] /  sales_mean
            sim.sales_d_sd_mean[t,1] = sim.sales_d_sd[t,1] /  sales_d_mean

            labor_mean= sum(real(measure.*(n+ (e.*F_mat./p_o_t.w )*(1-s.fcost_fgoods))))
            labor_sd=sqrt.(sum(real(measure.*( n+ (e.*F_mat./p_o_t.w)*(1-s.fcost_fgoods) .- labor_mean ).^2)))
            sim.labor_sd_mean[t,1] = labor_sd / labor_mean


           # Average productivity
            sim.avg_productivity1[t,1] = ((sum(real(measure.*(rt[t].z).^(σ-1))))/sum(measure) ).^(1/(σ - 1))
            sim.avg_productivity2[t,1] =  sum(real(measure.*sales.* rt[t].z)) / sum(real(measure .* sales))
            sim.avg_productivity3[t,1] =  sum(real(measure.*sales_d.* rt[t].z)) / sum(real(measure .* sales_d ))
        end


    ## Real Statistics
    #Real GDP, Real Exports and Real Domestic Sales

        sim.X_Laspeyres[t,1] = sum(real(measure.*p_o_t.ξt[1].*rt[1].pf.*rt[t].yf))
        sim.D_Laspeyres[t,1] = sum(real(measure.*rt[1].pd.*rt[t].yd))
        sim.Sales_Laspeyres[t,1] = sum(real(measure.*(rt[1].pd.*rt[t].yd+p_o_t.ξt[1].*rt[1].pf.*rt[t].yf)))
        sim.GDP_Laspeyres[t,1] = sum(real(measure.*(rt[1].pd.*rt[t].yd+p_o_t.ξt[1]*rt[1].pf.*rt[t].yf-mat*p_o_t.Pkt[1])))
        sim.Imports[t,1] = sim.PmYm[t,1]
        sim.Imports_C[t,1] = p_o_t.ξ*(1+τ_m_c)*Pm_c*ym_c
        sim.Imports_K[t,1] =  p_o_t.ξ*(1+τ_m_k)*Pm_k*ym_k
        sim.Imports_C_Laspeyres[t,1] = p_o_t.ξt[1].*(1+τ_m_c_v[1])*Pm_c*ym_c
        sim.Imports_K_Laspeyres[t,1] =  p_o_t.ξt[1].*(1+τ_m_k_v[1])*Pm_k*ym_k
        sim.Imports_Laspeyres[t,1] =sim.Imports_C_Laspeyres[t,1] +sim.Imports_K_Laspeyres[t,1]

     ## Solution-related statistics

        sim.a_min_share[t,1] = sum(real(measure[1,:]))
        sim.a_max_share[t,1] = sum(real(measure[end,:]))

    end
return sim,rt
end

############################# Transition_Vec2 (2 versions, first for using NLsolve [in-place function], second standard)#############################

####### For NLsolve #######
#This version of KLS3_transition_vec2 is in-place, that is, its output (market clearing conditions) is a modified value of an input (F, in this case). This is done because NLsolve needs this format to work. It is customary to add a ! at the end of in-place functions.

function KLS3_transition_vec2!(F,Guess,p,p_o_t,rt)
@unpack N, ω_h_c, σ, ω_h_k,rv,Pfv,Yfv,θ_v,β_v,δ_v,τ_v,τ_x_v,τ_m_c_v,τ_m_k_v,F_base,ω_m_k,ω_m_c,Pm_c,τ_m_c,τ_m_k,tariffsincome,Pm_k = p #See comment at the beginning of KL3_measure_trans

##Transform vector of guess into standard matrix form (used for LeastSquaresOptim solver)
if s.solver_LM==1
Guess1=copy(Guess)
Guess=hcat(Guess1[1:N-2],Guess1[N-1:2*N-4],Guess1[2*N-3:3*N-6],Guess1[3*N-5:4*N-8],Guess1[4*N-7:5*N-10])
end

	wguess = [p_o_t.wt[1] exp.(Guess[N-1:2N-4])' p_o_t.wt[N]]
	ξguess = [p_o_t.ξt[1] exp.(Guess[2N-3:3N-6])' p_o_t.ξt[N]]
	Pkguess = [p_o_t.Pkt[1] exp.(Guess[3N-5:4N-8])' p_o_t.Pkt[N]]

    if s.tariffsincome == 1
		Ycguess = [p_o_t.Yct[1] exp.(Guess[1:N-2])' p_o_t.Yct[N]]
        Ykguess = [p_o_t.Ykt[1] exp.(Guess[4*N-7:5*N-10])' p_o_t.Ykt[N]]

        ϕhguess= (ω_h_c.^σ).*Ycguess + ((ω_h_k.*Pkguess).^σ).*Ykguess

    else
		ϕhguess = [p_o_t.ϕht[1] exp.(Guess[1:N-2])' p_o_t.ϕht[N]]

    end

    for t=N-1:-1:2

        ## Shocks

		# shock to interest rate
        p_o_t=merge(p_o_t,(r = p.rv[t],
        r_lag = p.rv[t-1],

        # shock to foreign price
        Pf_lag = p.Pfv[t-1],
        Pf = p.Pfv[t],

        # shock to foreign demand
        Yf = p.Yfv[t],

        # shock to collateral constraint
        θ = p.θ_v[t],

        # shock to discount factor
        β = p.β_v[t],

        # shock to depreciation rate
        δ = p.δ_v[t],

        # shock to iceberg costs
        τ = p.τ_v[t],

        # shock to tariffs
        τ_x = p.τ_x_v[t],
        τ_m_c = p.τ_m_c_v[t],
        τ_m_k = p.τ_m_k_v[t],

        ## GE prices ##

        ϕ_h = ϕhguess[t],
        ϕ_h_lag = ϕhguess[t-1],

        w = wguess[t],
        w_lag = wguess[t-1],

        ξ = ξguess[t],
        ξ_lag = ξguess[t-1],

        Pk = Pkguess[t],
        Pk_lag = Pkguess[t-1]))


		## Fixed costs
        if  s.fcost_fgoods==0 # If in units of labor
            p_o_t=merge(p_o_t,(F = p_o_t.w*F_base,))
         else
            p_o_t=merge(p_o_t,(F = F_base,))
        end

      ## Tariff income (initialized to 0) ##
        if s.tariffsincome == 1

            p_o_t=merge(p_o_t,(Yk = Ykguess[t],
            Yc = Ycguess[t]))
            #Yc=p_o_t.Yc
            #Yc =  (ϕhguess[t] - (ω_m_k*Pkguess[t])^σ*Ykguess[t])./(ω_m_c^σ)
            ym_c = p_o_t.Yc*(p_o_t.ξ*Pm_c*(1+τ_m_c)/ω_m_c)^(-σ)
            ym_k = p_o_t.Yk*(Pkguess[t]^σ)*(p_o_t.ξ*Pm_k*(1+τ_m_k)/ω_m_k)^(-σ)
            if s.tariffsincome_PE==0
                p_o_t=merge(p_o_t,(tariffsincome = τ_m_c*p_o_t.ξ*Pm_c*ym_c + τ_m_k*p_o_t.ξ*Pm_k*ym_k,))
			else
				p_o_t=merge(p_o_t,(tariffsincome = tariffsincome,))
            end
            p_o_t.tariffsincomet[t]=p_o_t.tariffsincome
        end

        ## Solve static and dynamic problems
        # Period 2
		if t==2
			p_o_temp = KLS3_staticproblem_period2_altTiming(p_o_t,p)
		else
			p_o_temp = KLS3_staticproblem(p_o_t,p)
		end

        #Fixed costs
		p_o_temp=merge(p_o_temp,(S_mat = p_o_temp.e*p_o_t.F,
        F_mat = p_o_temp.e*p_o_t.F,
        profits=p_o_t.w .+ p_o_t.tariffsincome .+ p_o_temp.e.*p_o_temp.π_x .+ (1 .-p_o_temp.e).*p_o_temp.π_nx))
        rt[t]=p_o_temp

    end

	rt = KLS3_dynamicproblem_trans_vec_t(p_o_t,p,rt)


	if s.flag_simulate == 0
		sim_fun, rt = KLS3_simulate_trans(p,p_o_t,rt,Guess)
	end



    # Market clearing conditions

    if s.tariffsincome==0
        for j=1:N-2
        F[j]=sim_fun.mc_n[j+1]
        end
        for j=N-1:2*N-4
        F[j]=sim_fun.mc_y[j+1-N+2]
        end
        for j=2*N-3:3*N-6
        F[j]=sim_fun.mc_k[2:N-1]
        end
        for j=3*N-5:4*N-8
        F[j]=sim_fun.mc_y_belief[2:N-1]
        end
    elseif s.tariffsincome==1
        for j=1:N-2
        F[j]=sim_fun.mc_n[j+1]
        end
        for j=N-1:2*N-4
        F[j]=sim_fun.mc_y[j+3-N]
        end
        for j=2*N-3:3*N-6
        F[j]=sim_fun.mc_k[j+5-2*N]
        end
        for j=3*N-5:4*N-8
        F[j]=sim_fun.mc_y_belief[j+7-3*N]
        end
        for j=4*N-7:5*N-10
        F[j]=sim_fun.mc_Yk_belief[j+9-4*N]
        end
    end
    return F
end


####################### Standard #######################

function KLS3_transition_vec2(Guess,p,p_o_t,rt)
#See comments of the in-place version (KLS3_transition_vec2!)

@unpack N, ω_h_c, σ, ω_h_k,rv,Pfv,Yfv,θ_v,β_v,δ_v,τ_v,τ_x_v,τ_m_c_v,τ_m_k_v,F_base,ω_m_k,ω_m_c,Pm_c,τ_m_c,τ_m_k,tariffsincome,Pm_k=p

    wguess = [p_o_t.wt[1] exp.(Guess[N-1:2N-4])' p_o_t.wt[N]]
    ξguess = [p_o_t.ξt[1] exp.(Guess[2N-3:3N-6])' p_o_t.ξt[N]]
    Pkguess = [p_o_t.Pkt[1] exp.(Guess[3N-5:4N-8])' p_o_t.Pkt[N]]

    if s.tariffsincome == 1
        Ycguess = [p_o_t.Yct[1] exp.(Guess[1:N-2])' p_o_t.Yct[N]]
        Ykguess = [p_o_t.Ykt[1] exp.(Guess[4*N-7:5*N-10])' p_o_t.Ykt[N]]

        ϕhguess= (ω_h_c.^σ).*Ycguess + ((ω_h_k.*Pkguess).^σ).*Ykguess


    else
        ϕhguess = [p_o_t.ϕht[1] exp.(Guess[1:N-2])' p_o_t.ϕht[N]]

    end

    for t=N-1:-1:2

        ## Shocks

        #global p_o_t

        # shock to interest rate
        p_o_t=merge(p_o_t,(r = p.rv[t],
        r_lag = p.rv[t-1],

        # shock to foreign price
        Pf_lag = p.Pfv[t-1],
        Pf = p.Pfv[t],

        # shock to foreign demand
        Yf = p.Yfv[t],

        # shock to collateral constraint
        θ = p.θ_v[t],

        # shock to discount factor
        β = p.β_v[t],

        # shock to depreciation rate
        δ = p.δ_v[t],

        # shock to iceberg costs
        τ = p.τ_v[t],

        # shock to tariffs
        τ_x = p.τ_x_v[t],
        τ_m_c = p.τ_m_c_v[t],
        τ_m_k = p.τ_m_k_v[t],

        ## GE prices ##

        ϕ_h = ϕhguess[t],
        ϕ_h_lag = ϕhguess[t-1],

        w = wguess[t],
        w_lag = wguess[t-1],

        ξ = ξguess[t],
        ξ_lag = ξguess[t-1],

        Pk = Pkguess[t],
        Pk_lag = Pkguess[t-1]))

        ## Fixed costs $$
        if  s.fcost_fgoods==0 # If in units of labor
            p_o_t=merge(p_o_t,(F = p_o_t.w*F_base,))
         else
            p_o_t=merge(p_o_t,(F = F_base,))
        end

      ## Tariff income (initialized to 0) ##
        if s.tariffsincome == 1

            p_o_t=merge(p_o_t,(Yk = Ykguess[t],
            Yc = Ycguess[t]))
            #Yc=p_o_t.Yc
            #Yc =  (ϕhguess[t] - (ω_m_k*Pkguess[t])^σ*Ykguess[t])./(ω_m_c^σ)
            ym_c = p_o_t.Yc*(p_o_t.ξ*Pm_c*(1+τ_m_c)/ω_m_c)^(-σ)
            ym_k = p_o_t.Yk*(Pkguess[t]^σ)*(p_o_t.ξ*Pm_k*(1+τ_m_k)/ω_m_k)^(-σ)
            if s.tariffsincome_PE==0
                p_o_t=merge(p_o_t,(tariffsincome = τ_m_c*p_o_t.ξ*Pm_c*ym_c + τ_m_k*p_o_t.ξ*Pm_k*ym_k,))
			else
				p_o_t=merge(p_o_t,(tariffsincome = tariffsincome,))
            end
            p_o_t.tariffsincomet[t]=p_o_t.tariffsincome
        end

        ## Solve static and dynamic problems ##

        # Period 2
        if t==2
            p_o_temp = KLS3_staticproblem_period2_altTiming(p_o_t,p)
        else
            p_o_temp = KLS3_staticproblem(p_o_t,p)
        end

        #Fixed costs
        p_o_temp=merge(p_o_temp,(S_mat = p_o_temp.e*p_o_t.F,
        F_mat = p_o_temp.e*p_o_t.F,
        profits=p_o_t.w .+ p_o_t.tariffsincome .+ p_o_temp.e.*p_o_temp.π_x .+ (1 .-p_o_temp.e).*p_o_temp.π_nx))
        rt[t]=p_o_temp
    end

    rt = KLS3_dynamicproblem_trans_vec_t(p_o_t,p,rt)

    #ap_ind_tr_Julia=zeros(100,100,20)
    #for t=1:20
        #ap_ind_tr_Julia[:,:,t]=rt[t].ap_ind
    #end
    #file=matopen("ap_ind_tr_Julia.mat","w")
    #write(file,"ap_ind_tr_Julia",ap_ind_tr_Julia)
    #close(file)


    if s.flag_simulate == 0
        sim_fun, rt = KLS3_simulate_trans(p,p_o_t,rt,Guess)
    end

    # Market clearing conditions


    if s.tariffsincome==0
         mc = [sim_fun.mc_n[2:N-1]' sim_fun.mc_y[2:N-1]' sim_fun.mc_k[2:N-1]' sim_fun.mc_y_belief[2:N-1]']'
    elseif s.tariffsincome==1
        mc = [sim_fun.mc_n[2:N-1]' sim_fun.mc_y[2:N-1]' sim_fun.mc_k[2:N-1]' sim_fun.mc_y_belief[2:N-1]' sim_fun.mc_Yk_belief[2:N-1]']'
    end

    return mc, p_o_t, sim_fun, rt
end
