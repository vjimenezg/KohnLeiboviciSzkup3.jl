############################# measure_trans, dynamicproblem_trans_vec_t, dynamicproblem_trans_vec, simulate_trans, transition_vec2  ########################

########################## Measure Trans #############################
function KLS3_measure_trans(rt,N,M0)

    numZ = length(rt{1,1}.z_grid);
    numA = length(rt{1,1}.a_grid);


    M = zeros(numA,numZ,N);
    M(:,:,1) = M0;
    M(:,:,2) = M0; % in the second period news arrives after states are determined
    Mnew = M0;

    P=rt{1,1}.z_P;
    for n=2:N-1
        apt_ind = rt{n}.ap_ind;
        Mold = Mnew;
        Mnew = zeros(numA,numZ);

        for i=1:numA %Old asset state
            for j=1:numZ %Old productivity state

                Mnew(apt_ind(i,j),:) = Mnew(apt_ind(i,j),:) + Mold(i,j)*P(j,:);

            end
        end

        M(:,:,n+1) = Mnew;

    end
return M

end

############## Dynamic Problem Transition Vectorized t ###############

function KLS3_dynamicproblem_trans_vec_t(m,s,r,rt)


### Dynamic problem: value functions and policy functions

### Shocks

#Initialize solution objects
    v_p=rt{s.N,1}.v

    v_new = similar(v_p)
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
        profits = rt{t,1}.profits
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

        v_pt=v_new'
        #Store output from value function iteration
        rt{t,1}.v = v_new
        rt{t,1}.ap = ap
        rt{t,1}.ap_ind= ap_ind
        rt{t,1}.c = profits + rt{t,1}.a_grid_mat.*(1+m.r) -rt{t,1}.ap

    ## Store output to be used in simulation

        rt{t,1}.pd = (1-rt{t,1}.e).*rt{t,1}.pd_nx + rt{t,1}.e.*rt{t,1}.pd_x
        rt{t,1}.yd = (1-rt{t,1}.e).*rt{t,1}.yd_nx + rt{t,1}.e.*rt{t,1}.yd_x
        rt{t,1}.pf = rt{t,1}.e.*rt{t,1}.pf_x
        rt{t,1}.yf =rt{t,1}.e.*rt{t,1}.yf_x
        rt{t,1}.k = (1-rt{t,1}.e).*rt{t,1}.k_nx + rt{t,1}.e.*rt{t,1}.k_x
        rt{t,1}.n_x = rt{t,1}.e.*(rt{t,1}.n_x)
        rt{t,1}.n_nx = (1-rt{t,1}.e).*rt{t,1}.n_nx
        rt{t,1}.n = (1-rt{t,1}.e).*rt{t,1}.n_nx + rt{t,1}.e.*rt{t,1}.n_x
        rt{t,1}.m = (1-rt{t,1}.e).*rt{t,1}.m_nx + rt{t,1}.e.*rt{t,1}.m_x

        rt{t,1}.pd_nx = []; rt{t,1}.pd_x = []; rt{t,1}.yd_nx=[];  rt{t,1}.yd_x=[]; rt{t,1}.pf_x=[]; rt{t,1}.yf_x=[];
        rt{t,1}.k_nx =[]; rt{t,1}.k_x = []; rt{t,1}.n_nx=[]; rt{t,1}.n_x=[]; rt{t,1}.m_nx=[]; rt{t,1}.m_x=[];

        rt{t,1}.π = (1-rt{t,1}.e).*rt{t,1}.π_nx + rt{t,1}.e.*rt{t,1}.π_x
        rt{t,1}.a = rt{t,1}.a_grid_mat


    end

#     #Store output from value function iteration
#     r.v = v_new
#     r.ap = r.a_grid[ap_ind]
#     r.ap_ind= ap_ind
#     r.c = profits + r.a_grid_mat.*(1+m.r) -r.ap
#
# ## Store output to be used in simulation
#
#     r.e= r.e    #Should we multiply by sec_choice here and below?
#     r.pd = (1-r.e).*r.pd_nx + r.e.*r.pd_x
#     r.yd = (1-r.e).*r.yd_nx + r.e.*r.yd_x
#     r.pf = r.e.*r.pf_x
#     r.yf =r.e.*r.yf_x
#     r.k = (1-r.e).*r.k_nx + r.e.*r.k_x
#     r.n_x = r.e.*(r.n_x)
#     r.n_nx = (1-r.e).*r.n_nx
#     r.n = (1-r.e).*r.n_nx + r.e.*r.n_x
#     r.m = (1-r.e).*r.m_nx + r.e.*r.m_x
#
#     r.pd_nx = []; r.pd_x = []; r.yd_nx=[];  r.yd_x=[]; r.pf_x=[]; r.yf_x=[];
#     r.k_nx =[]; r.k_x = []; r.n_nx=[]; r.n_x=[]; r.m_nx=[]; r.m_x=[];
#
#     r.π = (1-r.e).*r.π_nx + r.e.*r.π_x
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
    profits = m.w + m.tariffsincome + r.e.*r.π_x + (1-r.e).*r.π_nx

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
    r= merge((v=v_new,
    ap = ap,
    ap_ind=ap_ind,
    c = profits .+ r.a_grid_mat.*(1+m.r) .- r.ap,

### Store output to be used in simulation

    e  =  r.e,    #Should we multiply by sec_choice here and below?
    pd = (1-r.e).*r.pd_nx + r.e.*r.pd_x,
    yd = (1-r.e).*r.yd_nx + r.e.*r.yd_x,
    pf = r.e.*r.pf_x,
    yf =r.e.*r.yf_x,
    k = (1-r.e).*r.k_nx + r.e.*r.k_x,
    n_x = r.e.*(r.n_x),
    n_nx = (1-r.e).*r.n_nx,
    n = (1-r.e).*r.n_nx + r.e.*r.n_x,
    m = (1-r.e).*r.m_nx + r.e.*r.m_x),r)

    r.pd_nx = []; r.pd_x = []; r.yd_nx=[];  r.yd_x=[]; r.pf_x=[]; r.yf_x=[];
    r.k_nx =[]; r.k_x = []; r.n_nx=[]; r.n_x=[]; r.m_nx=[]; r.m_x=[];

    r=merge((π = (1-r.e).*r.π_nx + r.e.*r.π_x,
    a = r.a_grid_mat,

    #Fixed costs
    S_mat = r.e*m.F,
    F_mat = r.e*m.F),r)

return r
end
############## Simulate Transition ###############


function KLS3_simulate_trans(m,s,rt,Guess)


    wguess = [m.wt[1] exp.(Guess[s.N-1:2*s.N-4]) m.wt[s.N]]
    ξguess = [m.ξt[1] exp.(Guess[2*s.N-3:3*s.N-6]) m.ξt[s.N]]
    Pkguess = [m.Pkt[1] exp.(Guess[3*s.N-5:4*s.N-8]) m.Pkt[s.N]]

    if s.tariffsincome == 1
        Ycguess = [m.Yct[1] exp.(Guess[1:s.N-2]) m.Yct[s.N]]
        Ykguess = [m.Ykt[1] exp.(Guess[4*s.N-7:5*s.N-10]) m.Ykt[s.N]]
        ϕhguess= (m.ω_h_c.^m.σ).*Ycguess + ((m.ω_h_k.*Pkguess).^m.σ).*Ykguess

    else
        Φhguess = [m.ϕht[1] exp.(Guess[1:s.N-2]) m.ϕht[s.N]]

    end

    N = s.N

    #Measure of firms
    M0 = rt{1}.measure
    sim=(measure = KLS3_measure_trans(rt,N,M0),)
         ### FROM HERE ONWARDS PENDING, WORKING ON KLS3_measure_trans ###
    for t=2:N
    ## Shocks

        m.r = m.rv(t); # shock to interest rate
        m.r_lag = m.rv(t-1);
        m.Pf_lag = m.Pfv(t-1); # shock to foreign price
        m.Pf = m.Pfv(t);
        m.Yf = m.Yfv(t); # shock to foreign demand
        m.theta = m.theta_v(t); # shock to collateral constraint
        m.beta = m.beta_v(t);    # shock to discount factor
        m.delta = m.delta_v(t); # shock to depreciation rate
        m.tau = m.tau_v(t); # shock to iceberg costs
        m.tau_x = m.tau_x_v(t); # shocks to tariffs
        m.tau_m_c = m.tau_m_c_v(t); # shocks to consumption tariffs
        m.tau_m_k = m.tau_m_k_v(t); # shocks to capital tariffs


        ## GE prices
        m.Phi_h = Phihguess(t);
        m.w = wguess(t);
        m.Xi = Xiguess(t);
        m.Pk = Pkguess(t);
        m.Pk_lag = Pkguess(t-1);

        if s.tariffsincome==1
            m.Yk = Ykguess(t);
            m.Yc = Ycguess(t);
            m.tariffsincome=m.tariffsincomet(t);
        end

        ## Useful objects

        n = rt{t,1}.n;
        mat = rt{t,1}.m;  # intermediates
        measure = sim.measure(:,:,t);  # measure
        pd = rt{t,1}.pd;
        yd = rt{t,1}.yd;
        pf = rt{t,1}.pf;
        yf = rt{t,1}.yf;
        e = rt{t,1}.e ;
        k = rt{t,1}.k;
        if t<N
            kp = rt{t+1,1}.k;
        end
        F_mat = rt{t,1}.F_mat;
        c = rt{t,1}.c;

        # Saving prices
        sim.w(t,1)=m.w;
        sim.Xi(t,1)=m.Xi;
        sim.Pk(t,1) = m.Pk;
        sim.Phi_h(t,1) = m.Phi_h;

        # Exporters and non-exporters
        sim.share_x(t,1) = sum(sum(measure.*e)); # share of exporters
        sim.share_d(t,1) = sum(sum(measure.*(1-e)));  #Mass of producers who are non-exporters

        ## Capital good
        sim.K(t,1) = sum(sum(measure.*k));
        if t>1 && t<N
            measurep =  sim.measure(:,:,t+1);
            sim.inv_agg(t,1) = sum(sum(measurep.*kp)) - (1-m.delta_v(t))*sum(sum(measure.*k));
        elseif t==N #For last period (t=N)
            # We assume that we reached steady state when t=N
            sim.inv_agg(t,1) = m.delta*sim.K(t,1);
        end

        sim.I(t,1) = max(sim.inv_agg(t,1),0); #Demand for new investment goods
        sim.M(t,1) = sum(sum(measure.*mat));

        yd_k = ((pd/(m.Pk*m.omega_h_k)).^(-m.sigma)) * (sim.I(t,1) + sim.M(t,1));
        ym_k = ((m.Xi*m.Pm_k*(1+m.tau_m_k)/(m.Pk*m.omega_m_k)).^(-m.sigma)) * (sim.I(t,1) + sim.M(t,1));

        sim.Yk(t,1) = ( sum(sum(measure.*m.omega_h_k.*(yd_k.^((m.sigma-1)/m.sigma)) )) + m.omega_m_k*(ym_k.^((m.sigma-1)/m.sigma)) )^(m.sigma/(m.sigma-1));
        sim.Pk(t,1) = ( sum(sum(measure.*(m.omega_h_k^m.sigma).*(pd.^(1-m.sigma)) )) + (m.omega_m_k^m.sigma)*((m.Xi*m.Pm_k*(1+m.tau_m_k)).^(1-m.sigma)) )^(1/(1-m.sigma));

        ## Consumption good

        sim.C(t,1) = sum(sum(measure.*c));

        #yd_c = ((pd/m.omega_h_c).^(-m.sigma)) * sim.C(t,1);
        yd_c = max(yd - yd_k,0.00001);
        ym_c = (m.Xi*m.Pm_c*(1+m.tau_m_c)/(m.omega_m_c)).^(-m.sigma) * sim.C(t,1);
        sim.Yc(t,1) = ( sum(sum(measure.*m.omega_h_c.*(yd_c.^((m.sigma-1)/m.sigma)))) + m.omega_m_c*(ym_c.^((m.sigma-1)/m.sigma)) )^(m.sigma/(m.sigma-1));


        ## Phi_h

        sim.Phi_h(t,1) = (m.omega_h_c^m.sigma)*sim.Yc(t,1) + ((m.omega_h_k*m.Pk)^m.sigma)*sim.Yk(t,1);


        ## Compute price and quantity indexes

        #Domestic sales
         sim.PdYd(t,1) =sum( measure(pd>0).*(pd(pd>0).*yd(pd>0)));

        #Exports
         sim.PxYx(t,1) = m.Xi*sum(measure(pf>0).*(pf(pf>0).*yf(pf>0)));
         sim.PxYx_USD(t,1) = sim.PxYx(t,1)/m.Xi;

         #Imports
        sim.PmYm(t,1) = m.Xi*(1+m.tau_m_c)*m.Pm_c*ym_c + m.Xi*(1+m.tau_m_k)*m.Pm_k*ym_k;

        # Tariffs Income (T)
        sim.tariffsincome(t,1)=m.Xi*m.tau_m_c*m.Pm_c*ym_c + m.Xi*m.tau_m_k*m.Pm_k*ym_k;
        if s.tariffsincome_PE==1
            sim.tariffsincome(t,1) = m.tariffsincome;
        end
        ## Compute aggregate variables needed to evaluate market clearing conditions

        #Labor and capital
        sim.FC(t,1) =  sum(sum(measure.*e.*F_mat));
        sim.N(t,1) = sum(sum(measure.*n))    +((sim.FC(t,1))/m.w)*(1-s.fcost_fgoods);

        ## Market clearing conditions

        #Labor
        sim.n_supply(t,1) = 1;
        sim.n_demand(t,1) = sim.N(t,1);
        #sim.mc_n(t,1) = ((1+sim.n_demand(t,1))/(1+sim.n_supply(t,1)))-1;
        sim.mc_n(t,1) = log(sim.n_demand(t,1)/sim.n_supply(t,1));

        #Assets
        #This market clearing condition is the result of reformulating the debt market clearing condition
        #In the original model, the sum of debt has to equal zero. Once the model is reformulated, this condition becomes the one below.
        sim.a_supply(t,1) = sum(sum(measure.*rt{t,1}.a));
        sim.a_demand(t,1) = sim.K(t,1);
        #sim.mc_a(t,1) = (1+sim.a_demand(t,1))/(1+sim.a_supply(t,1))-1;
        sim.mc_a(t,1) = log(sim.a_demand(t,1)/sim.a_supply(t,1));

        #Final goods
        sim.y_supply(t,1) = sim.Yc(t,1);
        sim.y_demand(t,1) = sim.C(t,1);
        sim.mc_y(t,1) = log(sim.y_demand(t,1)/sim.y_supply(t,1));

        # Capital goods
        sim.k_supply(t,1) = sim.Yk(t,1);
        sim.k_demand(t,1) = sim.I(t,1) + sim.M(t,1) + sim.FC(t,1)*s.fcost_fgoods;
        sim.mc_k(t,1) = log(sim.k_demand(t,1)/sim.k_supply(t,1));
#         sim.mc_k(t,1) =sim.k_demand(t,1)/sim.k_supply(t,1)-1;

        #Beliefs
        #To solve the entrepreneur's problem, they need to have a belief about the aggregate price and quantity indexes in the market
        #In equilibrium, these beliefs need to be consistent with the actual aggregate prices and quantities
        #sim.mc_y_belief(t,1) = log(sim.Yh(t,1)/m.Yh); sim.mc_y_belief = 10*log((1+sim.Ycpi_t)/(1+m.Yh));  sim.mc_y_belief(t,1) = ((1+sim.Yh(t,1))/(1+m.Yh))-1;


        if s.tariffsincome==1
            sim.mc_Yk_belief(t,1) = log(sim.Yk(t,1)/m.Yk);
            sim.mc_tariffs_belief(t,1) = log(sim.tariffsincome(t,1)/m.tariffsincome);
            sim.mc_y_belief(t,1) = ((1+sim.Yc(t,1))/(1+m.Yc))-1;
        else
            sim.mc_Yk_belief(t,1) = 0;
            sim.mc_tariffs_belief(t,1) = 0;
            sim.mc_y_belief(t,1) = log(sim.Phi_h(t,1)/Phihguess(t));
        end




        ## Main statistics

        sim.Sales(t,1) = sim.PdYd(t,1) + sim.PxYx(t,1); #sim.PxYx is already denominated in units of the domestic final good
        sim.GDP(t,1) = sim.Sales(t,1) -  sim.M(t,1)*m.Pk; #

        sim.X_GDP(t,1) = sim.PxYx(t,1)/sim.GDP(t,1);
        sim.D_GDP(t,1)  = sim.PdYd(t,1) /sim.GDP(t,1) ;

        sim.X_Sales(t,1) = sim.PxYx(t,1)/sim.Sales(t,1);
        sim.D_Sales(t,1) = sim.PdYd(t,1)/sim.Sales(t,1);

         sim.NX_GDP(t,1) = (sim.PxYx(t,1)-sim.PmYm(t,1))/sim.GDP(t,1);
         sim.NX_Sales(t,1) = (sim.PxYx(t,1)-sim.PmYm(t,1))/sim.Sales(t,1);

         sim.X_D(t,1) = sim.PxYx(t,1)/sim.PdYd(t,1) ;

        # Every credit statistic only for tradables
        rt{t,1}.d = (1+m.r).*(m.Pk_lag*k - rt{t,1}.a);

        sim.credit(t,1) = sum(sum(measure.*max(rt{t,1}.d,0))); # only entrepreneurs/firms (for workers r.d is negative)
        sim.credit_gdp(t,1) = sim.credit(t,1)/sim.GDP(t,1);
        sim.d_agg(t,1) = sum(sum(measure.*rt{t,1}.d)); # for both firms and workers
        sim.NFA_GDP(t,1) = -sim.d_agg(t,1)/sim.GDP(t,1);

        sim.k_wagebill(t,1) = m.Pk_lag*sim.K(t,1)/(m.w*sim.n_supply(t,1));

        # Sales
        sales_x =  m.Xi*rt{t,1}.pf.*rt{t,1}.yf; #r.pf is denominated in foreign currency, so we adjust it
        sales_d = rt{t,1}.pd.*rt{t,1}.yd;
        sales = sales_d+sales_x;
        x_share = sales_x./sales;

        sim.x_share_av(t,1) =  sum(nansum(measure.*e.* x_share )) / sum(sum(measure.*e));
        sim.x_share_av(isnan(sim.x_share_av(t,1)))=1;

        if s.extra_results==1
            ln_sales = log(sales);
            ln_sales_d =log(sales_d);
            ln_sales_ind = isnan(ln_sales)==0;
            ln_sales_d_ind = isnan(ln_sales_d)==0;

            ln_sales_mean = sum(sum(measure(ln_sales_ind).*ln_sales(ln_sales_ind) ));
            sim.ln_sales_sd(t,1) = sqrt( sum(sum(measure(ln_sales_ind) .* (ln_sales(ln_sales_ind) - ln_sales_mean).^2 )) );

            ln_sales_d_mean = sum(sum(measure(ln_sales_d_ind) .* ln_sales_d(ln_sales_d_ind)));
            sim.ln_sales_d_sd(t,1) = sqrt( sum(sum(measure(ln_sales_d_ind) .* (ln_sales_d(ln_sales_d_ind) - ln_sales_d_mean).^2 )) );

            sales_mean = sum(sum(measure.*sales ));
            sim.sales_sd(t,1) = sqrt( sum(sum(measure .* (sales - sales_mean).^2 )) );

            sales_d_mean = sum(sum(measure .* sales_d)) ;
            sim.sales_d_sd(t,1) = sqrt( sum(sum(measure .* (sales_d - sales_d_mean).^2 )) );

            sim.sales_avg_nx(t,1) = sum(sum(measure.*(1-e).*sales)) / sum(sum(measure.*(1-e)));
            sim.labor_avg_nx(t,1) = sum(sum(measure.*(1-e).*(n ))) / sum(sum(measure.*(1-e)));
            sim.sales_d_avg_nx(t,1) = sim.sales_avg_nx(t,1);

            sim.sales_avg_x(t,1) = sum(sum(measure.*e.*sales))/ sum(sum(measure.*e));
            sim.labor_avg_x(t,1) = ( (sim.FC(t,1)/m.w)*(1-s.fcost_fgoods) + sum(sum(measure.*e.*(n )) ) ) / sum(sum(measure.*e));

        # These include fixed costs
            sim.labor_tot(t,1) = sim.N(t,1);  # total labor demand for production
            sim.labor_tot_x(t,1) = (sim.FC(t,1)/m.w)*(1-s.fcost_fgoods)  + sum(sum(measure.*e.*(n )) ) ;
            sim.labor_tot_nx(t,1) =  sum(sum(measure.*(1-e).*(n)) ) ;
            sim.labor_tot_w(t,1) = sum(sum( measure));

            sim.labor_x_share(t,1) = sim.labor_tot_x(t,1)/sim.labor_tot(t,1);
            sim.labor_nx_share(t,1) = sim.labor_tot_nx(t,1)/sim.labor_tot(t,1);

            sim.sales_d_avg_x(t,1) = sum(sum(measure.*e.*sales_d))/ sum(sum(measure.*e));

            sim.xpremium_sales(t,1) = sim.sales_avg_x(t,1)/sim.sales_avg_nx(t,1);
            sim.xpremium_labor(t,1) = sim.labor_avg_x(t,1)/sim.labor_avg_nx(t,1);

            sim.xpremium_sales_d(t,1) = sim.sales_d_avg_x(t,1)/sim.sales_d_avg_nx(t,1);

            sim.ln_sales_sd_mean(t,1) = sim.ln_sales_sd(t,1) /  ln_sales_mean;
            sim.ln_sales_d_sd_mean(t,1) = sim.ln_sales_d_sd(t,1) /  ln_sales_d_mean;
            sim.sales_sd_mean(t,1) = sim.sales_sd(t,1) /  sales_mean;
            sim.sales_d_sd_mean(t,1) = sim.sales_d_sd(t,1) /  sales_d_mean;

            labor_mean= sum(sum(measure.*(n+ (e.*F_mat./m.w )*(1-s.fcost_fgoods) ))) ;
            labor_sd=sqrt( sum(sum(measure.*( n+ (e.*F_mat./m.w)*(1-s.fcost_fgoods) - labor_mean ).^2))  );
            sim.labor_sd_mean(t,1) = labor_sd /  labor_mean;


           # Average productivity
            sim.avg_productivity1(t,1) = (((sum(sum(measure.*(rt{t,1}.z).^(m.sigma-1) )))/sum(sum(measure)) ).^(1/(m.sigma - 1)));
            sim.avg_productivity2(t,1) =  sum(sum(measure.*sales.* rt{t,1}.z  )) / sum(sum(measure .* sales ));
            sim.avg_productivity3(t,1) =  sum(sum(measure.*sales_d.* rt{t,1}.z  )) / sum(sum(measure .* sales_d ));
        end


    ## Real Statistics
    #Real GDP, Real Exports and Real Domestic Sales

        sim.X_Laspeyres(t,1) = sum(sum(measure.*m.Xit(1).*rt{1}.pf.*rt{t}.yf));
        sim.D_Laspeyres(t,1) = sum(sum(measure.*rt{1}.pd.*rt{t}.yd));
        sim.Sales_Laspeyres(t,1) = sum(sum(measure.*(rt{1}.pd.*rt{t}.yd+m.Xit(1).*rt{1}.pf.*rt{t}.yf)));
        sim.GDP_Laspeyres(t,1) = sum(sum(measure.*(rt{1}.pd.*rt{t}.yd+m.Xit(1)*rt{1}.pf.*rt{t}.yf-mat*m.Pkt(1))));
        sim.Imports(t,1) = sim.PmYm(t,1);
        sim.Imports_C(t,1) = m.Xi*(1+m.tau_m_c)*m.Pm_c*ym_c;
        sim.Imports_K(t,1) =  m.Xi*(1+m.tau_m_k)*m.Pm_k*ym_k;
        sim.Imports_C_Laspeyres(t,1) = m.Xit(1).*(1+m.tau_m_c_v(1))*m.Pm_c*ym_c;
        sim.Imports_K_Laspeyres(t,1) =  m.Xit(1).*(1+m.tau_m_k_v(1))*m.Pm_k*ym_k;
        sim.Imports_Laspeyres(t,1) =sim.Imports_C_Laspeyres(t,1) +sim.Imports_K_Laspeyres(t,1);

     ## Solution-related statistics

        sim.a_min_share(t,1) = sum(measure(1,:));
        sim.a_max_share(t,1) = sum(measure(end,:));

    end
return sim,rt

end

####################### Transition Vectorized2 #######################

function KLS3_transition_vec2(Guess,m,r,s,rt)


    wguess = [m.wt[1] exp.(Guess[s.N-1:2*s.N-4]) m.wt[s.N]]
    ξguess = [m.ξt[1] exp.(Guess[2*s.N-3:3*s.N-6]) m.ξt[s.N]]
    Pkguess = [m.Pkt[1] exp.(Guess[3*s.N-5:4*s.N-8]) m.Pkt[s.N]]
    r_initial = r

    if s.tariffsincome == 1
        Ycguess = [m.Yct[1] exp.(Guess[1:s.N-2]) m.Yct[s.N]]
        Ykguess = [m.Ykt[1] exp.(Guess[4*s.N-7:5*s.N-10]) m.Ykt[s.N]]

        ϕhguess= (m.ω_h_c.^m.σ).*Ycguess + ((m.ω_h_k.*Pkguess).^m.σ).*Ykguess


    else
        ϕhguess = [m.ϕht[1] exp.(Guess[1:s.N-2]) m.ϕht[s.N]]

    end

    for t=s.N-1:-1:2

        ## Shocks

        # shock to interest rate
        m=merge((r = m.rv[t],
        r_lag = m.rv[t-1],

        # shock to foreign price
        Pf_lag = m.Pfv[t-1],
        Pf = m.Pfv[t],

        # shock to foreign demand
        Yf = Yfv[t],

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
        Pk_lag = Pkguess[t-1]),m)

        #Fixed costs
        if  s.fcost_fgoods==0 # If in units of labor
            m=merge((m.F = m.w*m.F_base,),m)
         else
            m=merge((m.F = m.F_base,),m)
        end

      # Tariff income (initialized to 0)
        if s.tariffsincome == 1

            m=merge((Yk = Ykguess[t],
            Yc = Ycguess[t]),m)
            #Yc=m.Yc
            #Yc =  (ϕhguess[t] - (m.ω_m_k*Pkguess[t])^m.σ*Ykguess[t])./(m.ω_m_c^m.σ)
            ym_c = m.Yc*(m.ξ*m.Pm_c*(1+m.τ_m_c)/m.ω_m_c)^(-m.σ)
            ym_k = m.Yk*(Pkguess[t]^m.σ)*(m.ξ*m.Pm_k*(1+m.τ_m_k)/m.ω_m_k)^(-m.σ)
            if s.tariffsincome_PE==0
                m=merge((tariffsincome = m.τ_m_c*m.ξ*m.Pm_c*ym_c + m.τ_m_k*m.ξ*m.Pm_k*ym_k,),m)
            end
            tariffsincomet[t]=m.tariffsincome
            m=merge((tariffsincomet=tariffsincomet,),m)
        end

        ## Solve static and dynamic problems

        # Period 2
        if t==2
            r_temp = KLS3_staticproblem_period2_altTiming(m,s,r_initial,rt)
        else
            r_temp = KLS3_staticproblem(m,s,r_initial)
        end

        #Fixed costs
        r_temp=merge((S_mat = r_temp.e*m.F,
        F_mat = r_temp.e*m.F,
        profits=m.w + m.tariffsincome + r_temp.e.*r_temp.π_x + (1-r_temp.e).*r_temp.π_nx),r_temp)
        rt{t,1}=r_temp

    end


    rt = KLS3_dynamicproblem_trans_vec_t(m,s,r,rt);

#         #Value function
#         vp=rt{t+1}.v,
#
#         #Dynamic problem and simulation (No sunk costs)
#         rt{t,1} = KLS3_dynamicproblem_trans_vec(m,s,r_temp,vp),
#
#     end

    if s.flag_simulate == 0



        sim_fun, rt = KLS3_simulate_trans(m,s,rt,Guess);

         ### FROM HERE ONWARDS PENDING, WORKING ON KLS3_simulate_trans ###

#         if s.tariffsincome==0
#             [sim_fun, rt] = KLS3_simulate_trans(m,s,rt,Guess);
#         else
#             [sim_fun, rt] = KLS3_simulate_trans(m,s,rt,Guess);
#         end

    elseif s.flag_simulate == 1


#          #save mat_temp_shocks;
#          load mat_temp_shocks;
#          s.N=5000000;
#          s.extra_results=2;

#         [sim_fun, rt] = KLS3_simulate_shock_trans(m,s,rt,Yguess,Xiguess,wguess,sim_0);

       #error('Simulation by generating random shocks not availabe for transition dynamics')

    end


    # A fix so that the code does not break
    sim_fun.mc_n(isnan(sim_fun.mc_n(2:s.N-1)==1)) = 10000;
    sim_fun.mc_y(isnan(sim_fun.mc_y(2:s.N-1))) = 10000;
    sim_fun.mc_y_belief(isnan(sim_fun.mc_y_belief(2:s.N-1))) = 10000;
    sim_fun.mc_y(isnan(sim_fun.mc_y(2:s.N-1))) = 10000;


    # Market clearing conditions


    if s.tariffsincome==0
         mc = [sim_fun.mc_n(2:s.N-1)' sim_fun.mc_y(2:s.N-1)' sim_fun.mc_k(2:s.N-1)' sim_fun.mc_y_belief(2:s.N-1)']';

    elseif s.tariffsincome==1

        #Fix so that code does not break
        sim_fun.mc_Yk_belief(isnan(sim_fun.mc_Yk_belief(2:s.N-1))) = 10000;

#         persistent count_trans
#         count_transif isempty(count_trans)
#             count_trans=1;
#         else
#             count_trans=count_trans+1;
#         end
#
#         if mod(count_trans,100)==0
#             disp(['count_trans: ' num2str(count_trans)]);
#         end

        mc = [sim_fun.mc_n(2:s.N-1)' sim_fun.mc_y(2:s.N-1)' sim_fun.mc_k(2:s.N-1)'...
              sim_fun.mc_y_belief(2:s.N-1)' sim_fun.mc_Yk_belief(2:s.N-1)']';

    end
    return mc, m, r, sim_fun, rt
end
