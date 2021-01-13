############################# measure_trans, dynamicproblem_trans_vec_t, dynamicproblem_trans_vec, simulate_trans, transition_vec2 (NLSolve), transition_vec2 (standard) ########################

########################## Measure Trans #############################
function KLS3_measure_trans(rt,N,M0)

    numZ = length(rt[1].r_0.z_grid)
    numA = length(rt[1].r_0.a_grid)


    M = zeros(numA,numZ,N)
    M[:,:,1] = M0
    M[:,:,2] = M0 # in the second period news arrives after states are determined
    Mnew = M0

    P=rt[1].r_0.z_P
    for n=2:N-1
        apt_ind = rt[n].ap_ind
        Mold = Mnew
        Mnew = zeros(numA,numZ)

        for i=1:numA #Old asset state
            for j=1:numZ #Old productivity state

                Mnew[apt_ind[i,j],:] = Mnew[apt_ind[i,j],:] .+ Mold[i,j].*P[j,:]

            end
        end

        M[:,:,n+1] = Mnew

    end
return M

end

############## Dynamic Problem Transition Vectorized t ###############

function KLS3_dynamicproblem_trans_vec_t(m,s,r,rt)


### Dynamic problem: value functions and policy functions

### Shocks

#Initialize solution objects
    v_p=rt[s.N].r_end.v

    v_new = fill(NaN,size(v_p))
    ap_ind_float = zeros(size(v_new))
    ap=ones(s.a_grid_size,s.z_grid_size)
    asset_income = r.a_grid'.*(1+m.r)
    exponent=1-m.γ
### Dynamic problem: value functions and policy functions
    P=r.z_P

    v_pt=v_p'
    a_grid_vec=r.a_grid_vec
    ones_vec=m.β *ones(s.a_grid_size,1)

    #Value function iteration algorithm


    for t=s.N-1:-1:2
        profits = rt[t].profits
        for j = 1:s.z_grid_size

                c =  profits[:,j] .+ asset_income .- a_grid_vec

                neg_c_indexes = c.<=0 #indices of a' for which consumption<0

                u = (c.^exponent)./exponent .+ neg_c_indexes*-1e50

                MM=ones_vec*P[j,:]'
                v,index_ap = findmax(u .+ MM*v_pt,dims=1)
                v_new[:,j] = v
                for (i,index) in enumerate(index_ap)
                ap_ind_float[i,j] = index[1]
                end
        end
        ap_ind=convert(Array{Int64},ap_ind_float)
        for j =1:s.z_grid_size
            ap[:,j]=r.a_grid[ap_ind[:,j]]
        end

        v_pt=v_new'
        #Store output from value function iteration
        rt[t]=merge(rt[t],(v = v_new,
        ap = ap,
        ap_ind= ap_ind))
        rt[t]=merge(rt[t],(c = profits + rt[t].a_grid_mat.*(1+m.r) .-rt[t].ap,

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
        a = rt[t].a_grid_mat))


    end

#     #Store output from value function iteration
#     r=merge(r,(v = v_new
#     ap = ap,
#     ap_ind= ap_ind,
#     c = profits .+ r.a_grid_mat.*(1+m.r) .- r.ap))
#
# ## Store output to be used in simulation
#
#     r=merge(r,(e= r.e,    #Should we multiply by sec_choice here and below?
#     pd = (1 .-r.e).*r.pd_nx + r.e.*r.pd_x,
#     yd = (1 .-r.e).*r.yd_nx + r.e.*r.yd_x,
#     pf = r.e.*r.pf_x,
#     yf =r.e.*r.yf_x,
#     k = (1 .-r.e).*r.k_nx + r.e.*r.k_x,
#     n_x = r.e.*(r.n_x),
#     n_nx = (1 .-r.e).*r.n_nx,
#     n = (1 .-r.e).*r.n_nx + r.e.*r.n_x,
#     m = (1 .-r.e).*r.m_nx + r.e.*r.m_x))
#
#     r=merge(r,(pd_nx = nothing, pd_x = nothing, yd_nx=nothing,  yd_x=nothing, pf_x=nothing, yf_x=nothing,
#     k_nx =nothing, k_x = nothing, n_nx=nothing, n_x=nothing, m_nx=nothing, m_x=nothing))
#
#     r.π_1 = (1 .-r.e).*r.π_nx + r.e.*r.π_x
#     r.a = r.a_grid_mat
#
#     #Fixed costs
#     r.S_mat = r.e*m.F
#     r.F_mat = r.e*m.F

return rt
end
############## Dynamic Problem Transition Vectorized ###############


function KLS3_dynamicproblem_trans_vec(m,s,r,v_p)

### Dynamic problem: value functions and policy functions

#Initialize solution objects
    M_one = ones(s.a_grid_size,s.z_grid_size)
    v_new = similar(v_p)
    ap_ind_float = zeros(size(v_new))
    ap=similar(M_one)
    asset_income = r.a_grid'.*(1+m.r)
    exponent=1-m.γ
### Dynamic problem: value functions and policy functions
    P=r.z_P
    v_pt=v_p'
    a_grid_vec=r.a_grid_vec
    ones_vec=m.β *ones(s.z_grid_size,1)

    #Value function iteration algorithm
    profits = m.w + m.tariffsincome + r.e.*r.π_x + (1 .-r.e).*r.π_nx

    for j = 1:s.z_grid_size

            c =  profits[:,j] .+ asset_income .- a_grid_vec

            neg_c_indexes = c.<=0 #indices of a' for which consumption<0

            u = (c.^exponent)./exponent .+ neg_c_indexes*-1e50

            MM=ones_vec*P[j,:]
            v,index_ap = findmax(u .+ MM*v_pt,dims=1)
            v_new[:,j] = v
            for (i,index) in enumerate(index_ap)
            ap_ind_float[i,j] = index[1]
            end
    end
    ap_ind=convert(Array{Int64},ap_ind_float)
    for j =1:s.z_grid_size
        ap[:,j]=r.a_grid[ap_ind[:,j]]
    end

    #Store output from value function iteration
    r= merge(r,(v=v_new,
    ap = ap,
    ap_ind=ap_ind,
    c = profits .+ r.a_grid_mat.*(1+m.r) .- r.ap,

### Store output to be used in simulation

    e  =  r.e,    #Should we multiply by sec_choice here and below?
    pd = (1 .-r.e).*r.pd_nx + r.e.*r.pd_x,
    yd = (1 .-r.e).*r.yd_nx + r.e.*r.yd_x,
    pf = r.e.*r.pf_x,
    yf =r.e.*r.yf_x,
    k = (1 .-r.e).*r.k_nx + r.e.*r.k_x,
    n_x = r.e.*(r.n_x),
    n_nx = (1 .-r.e).*r.n_nx,
    n = (1 .-r.e).*r.n_nx + r.e.*r.n_x,
    m = (1 .-r.e).*r.m_nx + r.e.*r.m_x))

    r=merge(r,(pd_nx = nothing, pd_x = nothing, yd_nx=nothing,  yd_x=nothing, pf_x=nothing, yf_x=nothing,
    k_nx =nothing, k_x = nothing, n_nx=nothing, n_x=nothing, m_nx=nothing, m_x=nothing))

    r=merge(r,(π_1 = (1 .-r.e).*r.π_nx + r.e.*r.π_x,
    a = r.a_grid_mat,

    #Fixed costs
    S_mat = r.e*m.F,
    F_mat = r.e*m.F))

return r
end
############## Simulate Transition ###############


function KLS3_simulate_trans(m,s,rt,Guess)


    wguess = [m.wt[1] exp.(Guess[s.N-1:2*s.N-4])' m.wt[s.N]]
    ξguess = [m.ξt[1] exp.(Guess[2*s.N-3:3*s.N-6])' m.ξt[s.N]]
    Pkguess = [m.Pkt[1] exp.(Guess[3*s.N-5:4*s.N-8])' m.Pkt[s.N]]

    if s.tariffsincome == 1
        Ycguess = [m.Yct[1] exp.(Guess[1:s.N-2])' m.Yct[s.N]]
        Ykguess = [m.Ykt[1] exp.(Guess[4*s.N-7:5*s.N-10])' m.Ykt[s.N]]
        ϕhguess= (m.ω_h_c.^m.σ).*Ycguess + ((m.ω_h_k.*Pkguess).^m.σ).*Ykguess

    else
        ϕhguess = [m.ϕht[1] exp.(Guess[1:s.N-2])' m.ϕht[s.N]]

    end

    N = s.N

    #Measure of firms
    M0 = rt[1].measure
    sim=(measure = KLS3_measure_trans(rt,N,M0),)

    w_sim=zeros(N,1)
    ξ_sim=zeros(N,1)
    Pk_sim=zeros(N,1)
    ϕ_h_sim=zeros(N,1)
    share_x_sim=zeros(N,1)
    share_d_sim=zeros(N,1)
    K_sim=zeros(N,1)
    inv_agg_sim=zeros(N,1)
    I_sim=zeros(N,1)
    M_sim=zeros(N,1)
    Yk_sim=zeros(N,1)
    C_sim=zeros(N,1)
    Yc_sim=zeros(N,1)
    PdYd_sim=zeros(N,1)
    PxYx_sim=zeros(N,1)
    PxYx_USD_sim=zeros(N,1)
    PmYm_sim=zeros(N,1)
    tariffsincome_sim=zeros(N,1)
    FC_sim=zeros(N,1)
    N_sim=zeros(N,1)
    n_supply_sim=zeros(N,1)
    n_demand_sim=zeros(N,1)
    mc_n_sim=zeros(N,1)
    a_supply_sim=zeros(N,1)
    a_demand_sim=zeros(N,1)
    mc_a_sim=zeros(N,1)
    y_supply_sim=zeros(N,1)
    y_demand_sim=zeros(N,1)
    mc_y_sim=zeros(N,1)
    k_supply_sim=zeros(N,1)
    k_demand_sim=zeros(N,1)
    mc_k_sim=zeros(N,1)
    mc_Yk_belief_sim=zeros(N,1)
    mc_tariffs_belief_sim=zeros(N,1)
    mc_y_belief_sim=zeros(N,1)
    Sales_sim=zeros(N,1)
    GDP_sim=zeros(N,1)
    X_GDP_sim=zeros(N,1)
    D_GDP_sim=zeros(N,1)
    X_Sales_sim=zeros(N,1)
    D_Sales_sim=zeros(N,1)
    NX_GDP_sim=zeros(N,1)
    NX_Sales_sim=zeros(N,1)
    X_D_sim=zeros(N,1)
    credit_sim=zeros(N,1)
    credit_gdp_sim=zeros(N,1)
    d_agg_sim=zeros(N,1)
    NFA_GDP_sim=zeros(N,1)
    k_wagebill_sim=zeros(N,1)
    x_share_av_sim=zeros(N,1)
    ln_sales_sd_sim=zeros(N,1)
    ln_sales_d_sd_sim=zeros(N,1)
    sales_sd_sim=zeros(N,1)
    sales_d_sd_sim=zeros(N,1)
    sales_avg_nx_sim=zeros(N,1)
    labor_avg_nx_sim=zeros(N,1)
    sales_d_avg_nx_sim=zeros(N,1)
    sales_avg_x_sim=zeros(N,1)
    labor_avg_x_sim=zeros(N,1)
    labor_tot_sim=zeros(N,1)
    labor_tot_x_sim=zeros(N,1)
    labor_tot_nx_sim=zeros(N,1)
    labor_tot_w_sim=zeros(N,1)
    labor_x_share_sim=zeros(N,1)
    labor_nx_share_sim=zeros(N,1)
    sales_d_avg_x_sim=zeros(N,1)
    xpremium_sales_sim=zeros(N,1)
    xpremium_labor_sim=zeros(N,1)
    xpremium_sales_d_sim=zeros(N,1)
    ln_sales_sd_mean_sim=zeros(N,1)
    ln_sales_d_sd_mean_sim=zeros(N,1)
    sales_sd_mean_sim=zeros(N,1)
    sales_d_sd_mean_sim=zeros(N,1)
    labor_sd_mean_sim=zeros(N,1)
    avg_productivity1_sim=zeros(N,1)
    avg_productivity2_sim=zeros(N,1)
    avg_productivity3_sim=zeros(N,1)
    X_Laspeyres_sim=zeros(N,1)
    D_Laspeyres_sim=zeros(N,1)
    Sales_Laspeyres_sim=zeros(N,1)
    GDP_Laspeyres_sim=zeros(N,1)
    Imports_sim=zeros(N,1)
    Imports_C_sim=zeros(N,1)
    Imports_K_sim=zeros(N,1)
    Imports_C_Laspeyres_sim=zeros(N,1)
    Imports_K_Laspeyres_sim=zeros(N,1)
    Imports_Laspeyres_sim=zeros(N,1)
    a_min_share_sim=zeros(N,1)
    a_max_share_sim=zeros(N,1)




    sim=merge(sim,(w=w_sim,ξ=ξ_sim,Pk=Pk_sim,ϕ_h=ϕ_h_sim,share_x=share_x_sim,share_d=share_d_sim,K=K_sim,inv_agg=inv_agg_sim,I=I_sim,M=M_sim,Yk=Yk_sim,C=C_sim,Yc=Yc_sim,PdYd=PdYd_sim,PxYx=PxYx_sim,PxYx_USD=PxYx_USD_sim,PmYm=PmYm_sim,tariffsincome=tariffsincome_sim,FC=FC_sim,N=N_sim,n_supply=n_supply_sim,n_demand=n_demand_sim,mc_n=mc_n_sim,a_supply=a_supply_sim,a_demand=a_demand_sim,mc_a=mc_a_sim,y_supply=y_supply_sim,y_demand=y_demand_sim,mc_y=mc_y_sim,k_supply=k_supply_sim,k_demand=k_demand_sim,mc_k=mc_k_sim,mc_Yk_belief=mc_Yk_belief_sim,mc_tariffs_belief=mc_tariffs_belief_sim,mc_y_belief=mc_y_belief_sim,Sales=Sales_sim,GDP=GDP_sim,X_GDP=X_GDP_sim,D_GDP=D_GDP_sim,X_Sales=X_Sales_sim,D_Sales=D_Sales_sim,NX_GDP=NX_GDP_sim,NX_Sales=NX_Sales_sim,X_D=X_D_sim,credit=credit_sim,credit_gdp=credit_gdp_sim,d_agg=d_agg_sim,NFA_GDP=NFA_GDP_sim,k_wagebill=k_wagebill_sim,x_share_av=x_share_av_sim,ln_sales_sd=ln_sales_sd_sim,ln_sales_d_sd=ln_sales_d_sd_sim,sales_sd=sales_sd_sim,sales_d_sd=sales_d_sd_sim,sales_avg_nx=sales_avg_nx_sim,labor_avg_nx=labor_avg_nx_sim,sales_d_avg_nx=sales_d_avg_nx_sim,sales_avg_x=sales_avg_x_sim,labor_avg_x=labor_avg_x_sim,labor_tot=labor_tot_sim,labor_tot_x=labor_tot_x_sim,labor_tot_nx=labor_tot_nx_sim,labor_tot_w=labor_tot_w_sim,labor_x_share=labor_x_share_sim,labor_nx_share=labor_nx_share_sim,sales_d_avg_x=sales_d_avg_x_sim,xpremium_sales=xpremium_sales_sim,xpremium_labor=xpremium_labor_sim,xpremium_sales_d=xpremium_sales_d_sim,ln_sales_sd_mean=ln_sales_sd_mean_sim,ln_sales_d_sd_mean=ln_sales_d_sd_mean_sim,sales_sd_mean=sales_sd_mean_sim,sales_d_sd_mean=sales_d_sd_mean_sim,labor_sd_mean=labor_sd_mean_sim,avg_productivity1=avg_productivity1_sim,avg_productivity2=avg_productivity2_sim,avg_productivity3=avg_productivity3_sim,X_Laspeyres=X_Laspeyres_sim,D_Laspeyres=D_Laspeyres_sim,Sales_Laspeyres=Sales_Laspeyres_sim,GDP_Laspeyres=GDP_Laspeyres_sim,Imports=Imports_sim,Imports_C=Imports_C_sim,Imports_K=Imports_K_sim,Imports_C_Laspeyres=Imports_C_Laspeyres_sim,Imports_K_Laspeyres=Imports_K_Laspeyres_sim,Imports_Laspeyres=Imports_Laspeyres_sim,a_min_share=a_min_share_sim,a_max_share=a_max_share_sim))


    for t=2:N
    ### Shocks

        m=merge(m,(r = m.rv[t], # shock to interest rate
        r_lag = m.rv[t-1],
        Pf_lag = m.Pfv[t-1], # shock to foreign price
        Pf = m.Pfv[t],
        Yf = m.Yfv[t], # shock to foreign demand
        θ = m.θ_v[t], # shock to collateral constraint
        β = m.β_v[t],    # shock to discount factor
        δ = m.δ_v[t], # shock to depreciation rate
        τ = m.τ_v[t], # shock to iceberg costs
        τ_x = m.τ_x_v[t], # shocks to tariffs
        τ_m_c = m.τ_m_c_v[t], # shocks to consumption tariffs
        τ_m_k = m.τ_m_k_v[t], # shocks to capital tariffs


        ## GE prices
        ϕ_h = ϕhguess[t],
        w = wguess[t],
        ξ = ξguess[t],
        Pk = Pkguess[t],
        Pk_lag = Pkguess[t-1]))

        if s.tariffsincome==1
            m=merge(m,(Yk = Ykguess[t],
            Yc = Ycguess[t],
            tariffsincome=m.tariffsincomet[t]))
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

        sim.w[t,1]=m.w
        sim.ξ[t,1]=m.ξ
        sim.Pk[t,1] = m.Pk
        sim.ϕ_h[t,1] = m.ϕ_h

        # Exporters and non-exporters
        sim.share_x[t,1] = sum(measure.*e) # share of exporters
        sim.share_d[t,1] = sum(measure.*(1 .-e))  #Mass of producers who are non-exporters

        ## Capital good
        sim.K[t,1] = sum(measure.*k)
        if t>1 && t<N
            measurep =  sim.measure[:,:,t+1]
            sim.inv_agg[t,1] = sum(measurep.*kp) - (1 .-m.δ_v[t])*sum(measure.*k)
        elseif t==N #For last period (t=N)
            # We assume that we reached steady state when t=N
            sim.inv_agg[t,1] = m.δ*sim.K[t,1]
        end

        sim.I[t,1] = max(sim.inv_agg[t,1],0) #Demand for new investment goods
        sim.M[t,1] = sum(measure.*mat)

        yd_k = ((pd/(m.Pk*m.ω_h_k)).^(-m.σ)) * (sim.I[t,1] + sim.M[t,1])
        ym_k = ((m.ξ*m.Pm_k*(1+m.τ_m_k)/(m.Pk*m.ω_m_k)).^(-m.σ)) * (sim.I[t,1] + sim.M[t,1])

        sim.Yk[t,1] = (sum(measure.*m.ω_h_k.*(yd_k.^((m.σ-1)/m.σ))) + m.ω_m_k*(ym_k.^((m.σ-1)/m.σ)))^(m.σ/(m.σ-1))
        sim.Pk[t,1] = (sum(measure.*(m.ω_h_k^m.σ).*(pd.^(1-m.σ))) + (m.ω_m_k^m.σ)*((m.ξ*m.Pm_k*(1+m.τ_m_k)).^(1-m.σ)) )^(1/(1-m.σ))

        ## Consumption good

        sim.C[t,1] = sum(measure.*c)

        #yd_c = ((pd/m.ω_h_c).^(-m.σ)) * sim.C[t,1]
        yd_c = max.(yd - yd_k,0.00001)
        ym_c = (m.ξ*m.Pm_c*(1+m.τ_m_c)/(m.ω_m_c)).^(-m.σ) * sim.C[t,1] # For some reason this has complex entries. Check dynamicproblem for errors in the computing of C
        sim.Yc[t,1] = (sum(measure.*m.ω_h_c.*(yd_c.^((m.σ-1)/m.σ))) + m.ω_m_c*(ym_c.^((m.σ-1)/m.σ)) )^(m.σ/(m.σ-1))


        ## ϕ_h

        sim.ϕ_h[t,1] = (m.ω_h_c^m.σ)*sim.Yc[t,1] + ((m.ω_h_k*m.Pk)^m.σ)*sim.Yk[t,1]


        ## Compute price and quantity indexes

        #Domestic sales
         sim.PdYd[t,1] =sum(measure[pd.>0].*(pd[pd.>0].*yd[pd.>0]))

        #Exports
         sim.PxYx[t,1] = m.ξ*sum(measure[pf.>0].*(pf[pf.>0].*yf[pf.>0]))

         sim.PxYx_USD[t,1] = sim.PxYx[t,1]/m.ξ

         #Imports
        sim.PmYm[t,1] = m.ξ*(1+m.τ_m_c)*m.Pm_c*ym_c + m.ξ*(1+m.τ_m_k)*m.Pm_k*ym_k

        # Tariffs Income (T)
        sim.tariffsincome[t,1]=m.ξ*m.τ_m_c*m.Pm_c*ym_c + m.ξ*m.τ_m_k*m.Pm_k*ym_k
        if s.tariffsincome_PE==1
            sim.tariffsincome[t,1] = m.tariffsincome
        end
        ## Compute aggregate variables needed to evaluate market clearing conditions

        #Labor and capital
        sim.FC[t,1] =  sum(measure.*e.*F_mat)
        sim.N[t,1] = sum(measure.*n)    +((sim.FC[t,1])/m.w)*(1 .-s.fcost_fgoods)

        ## Market clearing conditions

        #Labor
        sim.n_supply[t,1] = 1
        sim.n_demand[t,1] = sim.N[t,1]
        #sim.mc_n[t,1] = ((1+sim.n_demand[t,1])/(1+sim.n_supply[t,1]))-1
        sim.mc_n[t,1] = log.(sim.n_demand[t,1]/sim.n_supply[t,1])

        #Assets
        #This market clearing condition is the result of reformulating the debt market clearing condition
        #In the original model, the sum of debt has to equal zero. Once the model is reformulated, this condition becomes the one below.
        sim.a_supply[t,1] = sum(measure.*rt[t].a)
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
            sim.mc_Yk_belief[t,1] = log.(sim.Yk[t,1]/m.Yk)
            sim.mc_tariffs_belief[t,1] = log.(sim.tariffsincome[t,1]/m.tariffsincome)
            sim.mc_y_belief[t,1] = ((1+sim.Yc[t,1])/(1+m.Yc))-1
        else
            sim.mc_Yk_belief[t,1] = 0
            sim.mc_tariffs_belief[t,1] = 0
            sim.mc_y_belief[t,1] = log.(sim.ϕ_h[t,1]/ϕhguess[t])
        end




        ## Main statistics

        sim.Sales[t,1] = sim.PdYd[t,1] + sim.PxYx[t,1] #sim.PxYx is already denominated in units of the domestic final good
        sim.GDP[t,1] = sim.Sales[t,1] -  sim.M[t,1]*m.Pk #

        sim.X_GDP[t,1] = sim.PxYx[t,1]/sim.GDP[t,1]
        sim.D_GDP[t,1]  = sim.PdYd[t,1] /sim.GDP[t,1]

        sim.X_Sales[t,1] = sim.PxYx[t,1]/sim.Sales[t,1]
        sim.D_Sales[t,1] = sim.PdYd[t,1]/sim.Sales[t,1]

         sim.NX_GDP[t,1] = (sim.PxYx[t,1]-sim.PmYm[t,1])/sim.GDP[t,1]
         sim.NX_Sales[t,1] = (sim.PxYx[t,1]-sim.PmYm[t,1])/sim.Sales[t,1]

         sim.X_D[t,1] = sim.PxYx[t,1]/sim.PdYd[t,1]

        # Every credit statistic only for tradables
        rt[t].d = (1+m.r).*(m.Pk_lag*k - rt[t].a)

        sim.credit[t,1] = sum(measure.*max(rt[t].d,0)) # only entrepreneurs/firms (for workers r.d is negative)
        sim.credit_gdp[t,1] = sim.credit[t,1]/sim.GDP[t,1]
        sim.d_agg[t,1] = sum(measure.*rt[t].d) # for both firms and workers
        sim.NFA_GDP[t,1] = -sim.d_agg[t,1]/sim.GDP[t,1]

        sim.k_wagebill[t,1] = m.Pk_lag*sim.K[t,1]/(m.w*sim.n_supply[t,1])

        # Sales
        sales_x =  m.ξ*rt[t].pf.*rt[t].yf #r.pf is denominated in foreign currency, so we adjust it
        sales_d = rt[t].pd.*rt[t].yd
        sales = sales_d+sales_x
        x_share = sales_x./sales

        sim.x_share_av[t,1] =  NaNMath.sum(measure.*e.* x_share) / sum(measure.*e)
        sim.x_share_av[isnan.(sim.x_share_av[t,1])].=1

        if s.extra_results==1
            ln_sales = log.(sales)
            ln_sales_d =log.(sales_d)
            ln_sales_ind = isnan.(ln_sales).==0
            ln_sales_d_ind = isnan.(ln_sales_d).==0

            ln_sales_mean = sum(measure[ln_sales_ind].*ln_sales[ln_sales_ind])
            sim.ln_sales_sd[t,1] = sqrt.(sum(measure[ln_sales_ind] .* (ln_sales[ln_sales_ind] .- ln_sales_mean).^2 ))


            ln_sales_d_mean = sum(measure[ln_sales_d_ind] .* ln_sales_d[ln_sales_d_ind])
            sim.ln_sales_d_sd[t,1] = sqrt.( sum(measure[ln_sales_d_ind] .* (ln_sales_d[ln_sales_d_ind] .- ln_sales_d_mean).^2))

            sales_mean = sum(measure.*sales)
            sim.sales_sd[t,1] = sqrt.(sum(measure .* (sales .- sales_mean).^2))

            sales_d_mean = sum(measure .* sales_d)
            sim.sales_d_sd[t,1] = sqrt.(sum(measure .* (sales_d .- sales_d_mean).^2))

            sim.sales_avg_nx[t,1] = sum(measure.*(1 .-e).*sales) / sum(measure.*(1 .-e))
            sim.labor_avg_nx[t,1] = sum(measure.*(1 .-e).*(n)) / sum(measure.*(1 .-e))
            sim.sales_d_avg_nx[t,1] = sim.sales_avg_nx[t,1]

            sim.sales_avg_x[t,1] = sum(measure.*e.*sales)/ sum(measure.*e)
            sim.labor_avg_x[t,1] = ((sim.FC[t,1]/m.w)*(1-s.fcost_fgoods) + sum(measure.*e.*(n))) / sum(measure.*e)

        # These include fixed costs
            sim.labor_tot[t,1] = sim.N[t,1]  # total labor demand for production
            sim.labor_tot_x[t,1] = (sim.FC[t,1]/m.w)*(1-s.fcost_fgoods)  + sum(measure.*e.*(n))
            sim.labor_tot_nx[t,1] =  sum(measure.*(1 .-e).*(n))
            sim.labor_tot_w[t,1] = sum(measure)

            sim.labor_x_share[t,1] = sim.labor_tot_x[t,1]/sim.labor_tot[t,1]
            sim.labor_nx_share[t,1] = sim.labor_tot_nx[t,1]/sim.labor_tot[t,1]

            sim.sales_d_avg_x[t,1] = sum(measure.*e.*sales_d)/ sum(measure.*e)

            sim.xpremium_sales[t,1] = sim.sales_avg_x[t,1]/sim.sales_avg_nx[t,1]
            sim.xpremium_labor[t,1] = sim.labor_avg_x[t,1]/sim.labor_avg_nx[t,1]

            sim.xpremium_sales_d[t,1] = sim.sales_d_avg_x[t,1]/sim.sales_d_avg_nx[t,1]

            sim.ln_sales_sd_mean[t,1] = sim.ln_sales_sd[t,1] /  ln_sales_mean
            sim.ln_sales_d_sd_mean[t,1] = sim.ln_sales_d_sd[t,1] /  ln_sales_d_mean
            sim.sales_sd_mean[t,1] = sim.sales_sd[t,1] /  sales_mean
            sim.sales_d_sd_mean[t,1] = sim.sales_d_sd[t,1] /  sales_d_mean

            labor_mean= sum(measure.*(n+ (e.*F_mat./m.w )*(1-s.fcost_fgoods)))
            labor_sd=sqrt.(sum(measure.*( n+ (e.*F_mat./m.w)*(1-s.fcost_fgoods) .- labor_mean ).^2))
            sim.labor_sd_mean[t,1] = labor_sd /  labor_mean


           # Average productivity
            sim.avg_productivity1[t,1] = ((sum(measure.*(rt[t].z).^(m.sigma-1)))/sum(measure) ).^(1/(m.σ - 1))
            sim.avg_productivity2[t,1] =  sum(measure.*sales.* rt[t].z) / sum(measure .* sales)
            sim.avg_productivity3[t,1] =  sum(measure.*sales_d.* rt[t].z) / sum(measure .* sales_d )
        end


    ## Real Statistics
    #Real GDP, Real Exports and Real Domestic Sales

        sim.X_Laspeyres[t,1] = sum(measure.*m.ξt[1].*rt[1].r_0.pf.*rt[t].yf)
        sim.D_Laspeyres[t,1] = sum(measure.*rt[1].r_0.pd.*rt[t].yd)
        sim.Sales_Laspeyres[t,1] = sum(measure.*(rt[1].r_0.pd.*rt[t].yd+m.ξt[1].*rt[1].r_0.pf.*rt[t].yf))
        sim.GDP_Laspeyres[t,1] = sum(measure.*(rt[1].r_0.pd.*rt[t].yd+m.ξ_t[1]*rt[1].r_0.pf.*rt[t].yf-mat*m.Pkt[1]))
        sim.Imports[t,1] = sim.PmYm[t,1]
        sim.Imports_C[t,1] = m.ξ*(1+m.τ_m_c)*m.Pm_c*ym_c
        sim.Imports_K[t,1] =  m.ξ*(1+m.τ_m_k)*m.Pm_k*ym_k
        sim.Imports_C_Laspeyres[t,1] = m.ξt[1].*(1+m.τ_m_c_v(1))*m.Pm_c*ym_c
        sim.Imports_K_Laspeyres[t,1] =  m.ξt[1].*(1+m.τ_m_k_v(1))*m.Pm_k*ym_k
        sim.Imports_Laspeyres[t,1] =sim.Imports_C_Laspeyres[t,1] +sim.Imports_K_Laspeyres[t,1]

     ## Solution-related statistics

        sim.a_min_share[t,1] = sum(measure[1,:])
        sim.a_max_share[t,1] = sum(measure[end,:])

    end
return sim,rt
end

############################# Transition_Vec2 (2 versions, first for using NLsolve [in-place function], second standard)#############################
function KLS3_transition_vec2!(F,Guess,m,r,s,rt)


    wguess = [m.wt[1] exp.(Guess[s.N-1:2s.N-4])' m.wt[s.N]]
    ξguess = [m.ξt[1] exp.(Guess[2s.N-3:3s.N-6])' m.ξt[s.N]]
    Pkguess = [m.Pkt[1] exp.(Guess[3s.N-5:4s.N-8])' m.Pkt[s.N]]
    r_initial = r

    if s.tariffsincome == 1
        Ycguess = [m.Yct[1] exp.(Guess[1:s.N-2])' m.Yct[s.N]]
        Ykguess = [m.Ykt[1] exp.(Guess[4*s.N-7:5*s.N-10])' m.Ykt[s.N]]

        ϕhguess= (m.ω_h_c.^m.σ).*Ycguess + ((m.ω_h_k.*Pkguess).^m.σ).*Ykguess


    else
        ϕhguess = [m.ϕht[1] exp.(Guess[1:s.N-2])' m.ϕht[s.N]]

    end

    for t=s.N-1:-1:2

        ## Shocks

        # shock to interest rate
        m=merge(m,(r = m.rv[t],
        r_lag = m.rv[t-1],

        # shock to foreign price
        Pf_lag = m.Pfv[t-1],
        Pf = m.Pfv[t],

        # shock to foreign demand
        Yf = m.Yfv[t],

        # shock to collateral constraint
        θ = m.θ_v[t],

        # shock to discount factor
        β = m.β_v[t],

        # shock to depreciation rate
        δ = m.δ_v[t],

        # shock to iceberg costs
        τ = m.τ_v[t],

        # shock to tariffs
        τ_x = m.τ_x_v[t],
        τ_m_c = m.τ_m_c_v[t],
        τ_m_k = m.τ_m_k_v[t],


        ## GE prices

        ϕ_h = ϕhguess[t],
        ϕ_h_lag = ϕhguess[t-1],

        w = wguess[t],
        w_lag = wguess[t-1],

        ξ = ξguess[t],
        ξ_lag = ξguess[t-1],

        Pk = Pkguess[t],
        Pk_lag = Pkguess[t-1]))

        #Fixed costs
        if  s.fcost_fgoods==0 # If in units of labor
            m=merge(m,(F = m.w*m.F_base,))
         else
            m=merge(m,(F = m.F_base,))
        end

      # Tariff income (initialized to 0)

        if s.tariffsincome == 1
            m=merge(m,(Yk = Ykguess[t],
            Yc = Ycguess[t]))
            #Yc=m.Yc
            #Yc =  (ϕhguess[t] - (m.ω_m_k*Pkguess[t])^m.σ*Ykguess[t])./(m.ω_m_c^m.σ)
            ym_c = m.Yc*(m.ξ*m.Pm_c*(1+m.τ_m_c)/m.ω_m_c)^(-m.σ)
            ym_k = m.Yk*(Pkguess[t]^m.σ)*(m.ξ*m.Pm_k*(1+m.τ_m_k)/m.ω_m_k)^(-m.σ)
            if s.tariffsincome_PE==0
                m=merge(m,(tariffsincome = m.τ_m_c*m.ξ*m.Pm_c*ym_c + m.τ_m_k*m.ξ*m.Pm_k*ym_k,))
            end
            m.tariffsincomet[t]=m.tariffsincome
        end

        ## Solve static and dynamic problems

        # Period 2
        if t==2
            r_temp = KLS3_staticproblem_period2_altTiming(m,s,r_initial,rt)
        else
            r_temp = KLS3_staticproblem(m,s,r_initial)
        end

        #Fixed costs
        r_temp=merge(r_temp,(S_mat = r_temp.e*m.F,
        F_mat = r_temp.e*m.F,
        profits=m.w .+ m.tariffsincome .+ r_temp.e.*r_temp.π_x .+ (1 .-r_temp.e).*r_temp.π_nx))
        rt[t]=r_temp

    end

    rt = KLS3_dynamicproblem_trans_vec_t(m,s,r,rt)

    #         #Value function
    #         vp=rt[t+1].v,
    #
    #         #Dynamic problem and simulation (No sunk costs)
    #         rt[t] = KLS3_dynamicproblem_trans_vec(m,s,r_temp,vp),
    #
    #       end

    if s.flag_simulate == 0

        sim_fun, rt = KLS3_simulate_trans(m,s,rt,Guess)

    elseif s.flag_simulate == 1

        # THIS DOESN'T WORK (NO FUNCTION KLS3_simulate_shock_trans available, didn't translate the save/load either)
        #    #save mat_temp_shocks
        #         load mat_temp_shocks
        #         s.N=5000000
        #         s.extra_results=2

        #         sim_fun, rt = KLS3_simulate_shock_trans(m,s,rt,Yguess,ξguess,wguess,sim_0)

               #error('Simulation by generating random shocks not availabe for transition dynamics')
    end


    # A fix so that the code does not break
    sim_fun.mc_n[isnan.(sim_fun.mc_n)] .= 10000
    sim_fun.mc_y[isnan.(sim_fun.mc_y)] .= 10000
    sim_fun.mc_y_belief[isnan.(sim_fun.mc_y_belief)] .= 10000
    sim_fun.mc_y[isnan.(sim_fun.mc_y)] .= 10000


    # Market clearing conditions

    if s.tariffsincome==0
        for j=1:N-2
        F[j]=sim_fun.mc_n[j+1]
        end
        for j=N-1:2N-4
        F[j]=sim_fun.mc_y[j+1-N+2]
        end
        for j=2N-3:3N-6
        F[j]=sim_fun.mc_k[2:s.N-1]
        end
        for j=3N-5:4N-8
        F[j]=sim_fun.mc_y_belief[2:s.N-1]
        end
    elseif s.tariffsincome==1

        #Fix so that code does not break
        sim_fun.mc_Yk_belief[isnan.(sim_fun.mc_Yk_belief)] .= 10000
        for j=1:N-2
        F[j]=sim_fun.mc_n[j+1]
        end
        for j=N-1:2N-4
        F[j]=sim_fun.mc_y[j+3-N]
        end
        for j=2N-3:3N-6
        F[j]=sim_fun.mc_k[j+5-2N]
        end
        for j=3N-5:4N-8
        F[j]=sim_fun.mc_y_belief[j+7-3N]
        end
        for j=4N-7:5N-10
        F[j]=sim_fun.mc_Yk_belief[j+9-4N]
        end
    end
    return F
end


####################### Standard #######################

function KLS3_transition_vec2(Guess,m,r,s,rt)


    wguess = [m.wt[1] exp.(Guess[s.N-1:2s.N-4])' m.wt[s.N]]
    ξguess = [m.ξt[1] exp.(Guess[2s.N-3:3s.N-6])' m.ξt[s.N]]
    Pkguess = [m.Pkt[1] exp.(Guess[3s.N-5:4s.N-8])' m.Pkt[s.N]]
    r_initial = r

    if s.tariffsincome == 1
        Ycguess = [m.Yct[1] exp.(Guess[1:s.N-2])' m.Yct[s.N]]
        Ykguess = [m.Ykt[1] exp.(Guess[4*s.N-7:5*s.N-10])' m.Ykt[s.N]]

        ϕhguess= (m.ω_h_c.^m.σ).*Ycguess + ((m.ω_h_k.*Pkguess).^m.σ).*Ykguess


    else
        ϕhguess = [m.ϕht[1] exp.(Guess[1:s.N-2])' m.ϕht[s.N]]

    end

    for t=s.N-1:-1:2

        ## Shocks

        # shock to interest rate
        m=merge(m,(r = m.rv[t],
        r_lag = m.rv[t-1],

        # shock to foreign price
        Pf_lag = m.Pfv[t-1],
        Pf = m.Pfv[t],

        # shock to foreign demand
        Yf = m.Yfv[t],

        # shock to collateral constraint
        θ = m.θ_v[t],

        # shock to discount factor
        β = m.β_v[t],

        # shock to depreciation rate
        δ = m.δ_v[t],

        # shock to iceberg costs
        τ = m.τ_v[t],

        # shock to tariffs
        τ_x = m.τ_x_v[t],
        τ_m_c = m.τ_m_c_v[t],
        τ_m_k = m.τ_m_k_v[t],


        ## GE prices

        ϕ_h = ϕhguess[t],
        ϕ_h_lag = ϕhguess[t-1],

        w = wguess[t],
        w_lag = wguess[t-1],

        ξ = ξguess[t],
        ξ_lag = ξguess[t-1],

        Pk = Pkguess[t],
        Pk_lag = Pkguess[t-1]))

        #Fixed costs
        if  s.fcost_fgoods==0 # If in units of labor
            m=merge(m,(F = m.w*m.F_base,))
         else
            m=merge(m,(F = m.F_base,))
        end

      # Tariff income (initialized to 0)
        if s.tariffsincome == 1

            m=merge(m,(Yk = Ykguess[t],
            Yc = Ycguess[t]))
            #Yc=m.Yc
            #Yc =  (ϕhguess[t] - (m.ω_m_k*Pkguess[t])^m.σ*Ykguess[t])./(m.ω_m_c^m.σ)
            ym_c = m.Yc*(m.ξ*m.Pm_c*(1+m.τ_m_c)/m.ω_m_c)^(-m.σ)
            ym_k = m.Yk*(Pkguess[t]^m.σ)*(m.ξ*m.Pm_k*(1+m.τ_m_k)/m.ω_m_k)^(-m.σ)
            if s.tariffsincome_PE==0
                m=merge(m,(tariffsincome = m.τ_m_c*m.ξ*m.Pm_c*ym_c + m.τ_m_k*m.ξ*m.Pm_k*ym_k,))
            end
            m.tariffsincomet[t]=m.tariffsincome
        end

        ## Solve static and dynamic problems

        # Period 2
        if t==2
            r_temp = KLS3_staticproblem_period2_altTiming(m,s,r_initial,rt)
        else
            r_temp = KLS3_staticproblem(m,s,r_initial)
        end

        #Fixed costs
        r_temp=merge(r_temp,(S_mat = r_temp.e*m.F,
        F_mat = r_temp.e*m.F,
        profits=m.w .+ m.tariffsincome .+ r_temp.e.*r_temp.π_x .+ (1 .-r_temp.e).*r_temp.π_nx))
        rt[t]=r_temp

    end

    rt = KLS3_dynamicproblem_trans_vec_t(m,s,r,rt)

    #         #Value function
    #         vp=rt[t+1].v,
    #
    #         #Dynamic problem and simulation (No sunk costs)
    #         rt[t] = KLS3_dynamicproblem_trans_vec(m,s,r_temp,vp),
    #
    #     end

    if s.flag_simulate == 0

        sim_fun, rt = KLS3_simulate_trans(m,s,rt,Guess)

    elseif s.flag_simulate == 1

        # THIS DOESN'T WORK (NO FUNCTION KLS3_simulate_shock_trans available, didn't translate the save/load either)
        #    #save mat_temp_shocks
        #         load mat_temp_shocks
        #         s.N=5000000
        #         s.extra_results=2

        #         sim_fun, rt = KLS3_simulate_shock_trans(m,s,rt,Yguess,ξguess,wguess,sim_0)

               #error('Simulation by generating random shocks not availabe for transition dynamics')
    end


    # A fix so that the code does not break
    sim_fun.mc_n[isnan.(sim_fun.mc_n)] .= 10000
    sim_fun.mc_y[isnan.(sim_fun.mc_y)] .= 10000
    sim_fun.mc_y_belief[isnan.(sim_fun.mc_y_belief)] .= 10000
    sim_fun.mc_y[isnan.(sim_fun.mc_y)] .= 10000

    # Market clearing conditions


    if s.tariffsincome==0
         mc = [sim_fun.mc_n[2:s.N-1]' sim_fun.mc_y[2:s.N-1]' sim_fun.mc_k[2:s.N-1]' sim_fun.mc_y_belief[2:s.N-1]']'

    elseif s.tariffsincome==1

        #Fix so that code does not break
        sim_fun.mc_Yk_belief[isnan.(sim_fun.mc_Yk_belief)] .= 10000


        mc = [sim_fun.mc_n[2:s.N-1]' sim_fun.mc_y[2:s.N-1]' sim_fun.mc_k[2:s.N-1]' sim_fun.mc_y_belief[2:s.N-1]' sim_fun.mc_Yk_belief[2:s.N-1]']'

    end

    return mc, m, r, sim_fun, rt
end
