    ###############################################
    # Functions used in the computation of the Steady States
    # [Measure, Static Problem, Dynamic Problem, Simulation, GE Par (NLSolve), GE Par (Standard)]#####################################

############################# Measure ###########################

    function KLS3_measure(p_o,p,s,guessM)
    @unpack z_grid,a_grid,z_P=p # allows to use the parameters in p without having to addend "p.". Furthermore, it makes clear which parameters are used in each function.

    ## Computational objects
    numZ = length(z_grid)
    numA = length(a_grid)

    #Productivity
    P = z_P

    Mnew=guessM

    # Only valid with value function iteration
    ap_ind = p_o.ap_ind

    iter = 0
    diff_M = 1.0 #note that 1.0 =/= 1. The first is taken to be of type float64 while the second of type Int64. When diff_M is updated, it's gonna take a value that's not an exact integer, hence, it's better to declare it as float64 from the start (otherwise it changes type, something that makes code slower)

    while diff_M>s.eps_sm

        Mold = Mnew
        Mnew = zeros(numA,numZ)

        #for i=1:numA #Old asset state

        for j=1:numZ #Old productivity state
            #@inbounds begin
            PP=@view P[j,:]
            #aa_v=@view ap_ind[:,j]
                for i=1:numA #Old asset state
                #aa= ap_ind[i,j]
                #aa= aa_v[i]
                @views Mnew[ap_ind[i,j],:] = Mnew[ap_ind[i,j],:] + Mold[i,j]*PP #@views (and @view) lets P[j,:] reference the variable P instead of allocating a new array, P[j,:]. This speeds up the code and is thus used several times.
                end
            #end
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

function KLS3_staticproblem(p_o,p)

@unpack a_grid, z_grid_size,z,σ,r,δ,θ,Yf,τ,τ_x,α,α_m,F_base=p #See comment at the beginning of KL3_measure
### Useful objects

        #Assets
        p_o=merge(p_o,(a_grid_mat = a_grid'*ones(1,z_grid_size),

        #Productivity
        z_grid_mat = z))

        #ϕ_h

        ϕ_h = p_o.ϕ_h #Guessed value

        #Objects common for X and NX
        cap_gain = (1-p_o.Pk/p_o.Pk_lag)
        const_σ = ((σ-1)/σ)^σ
        p_o=merge(p_o,(rtilde_u = max.(((r+δ)+(1-δ)*cap_gain)*p_o.Pk_lag,0.00000001),))
        #note that without the max operator, rtilde_u takes negative values, making some variables below take complex values. This, in turn, makes the code break (because in Julia to deal with Complex numbers you gotta explicitly tell the compiler they are complex). As clearly those numbers can never be relevant, an easy fix is considering this maximum.
        if θ<1+r
            p_o=merge(p_o,(k_const = (1/p_o.Pk_lag)*((1+r)/(1+r-θ)).*p_o.a_grid_mat,))
        else
            p_o=merge(p_o,(k_const = ones(size(p_o.a_grid_mat)),))
        end

    ## Exporters
        ϕ_x = ϕ_h + (τ^(1-σ))*Yf*((p_o.ξ/(1+τ_x))^σ)

    #Unconstrained

        # marginal cost

        p_o=merge(p_o,(μ_u = (1 ./p_o.z_grid_mat) .* ((p_o.Pk/α_m).^α_m) .* ((p_o.w/((1-α)*(1-α_m))).^((1-α)*(1-α_m))) .* ((p_o.rtilde_u./(α*(1-α_m))).^(α*(1-α_m))),))


        p_o=merge(p_o,(k_x_u = (α*(1-α_m)./p_o.rtilde_u).*const_σ*ϕ_x .*(p_o.μ_u.^(1-σ)),
        n_x_u = ((1-α)*(1-α_m)/p_o.w)*const_σ*ϕ_x.*(p_o.μ_u.^(1-σ)),
        m_x_u = (α_m/p_o.Pk)*const_σ*ϕ_x.*(p_o.μ_u.^(1-σ)),

        yd_x_u = const_σ*ϕ_h*(p_o.μ_u.^(-σ)),
        yf_x_u = const_σ*(p_o.ξ^σ)*Yf*((1+τ_x)^(-σ))*(τ^(-σ)).*(p_o.μ_u.^(-σ)),

        pd_x_u = σ/(σ-1).*p_o.μ_u,
        pf_x_u = σ/(σ-1)*τ/p_o.ξ.*p_o.μ_u))

        #Solution
        p_o=merge(p_o,(k_x = p_o.k_x_u,
        n_x = p_o.n_x_u,
        m_x = p_o.m_x_u,
        yd_x = p_o.yd_x_u,
        yf_x = p_o.yf_x_u,
        pd_x = p_o.pd_x_u,
        pf_x = p_o.pf_x_u,
        const_x = BitArray(undef,size(p_o.k_x_u)))) #BitArray: Array of only 1's and 0's, uses less memory than standard Float64 arrays (hence, makes the code faster)

    #Constrained
        if θ<1+r

            p_o=merge(p_o,(rtilde_x_c = (((1 ./p_o.k_const) .* α * (1-α_m) * const_σ * ϕ_x) .* ((p_o.μ_u .* (p_o.rtilde_u).^(-α*(1-α_m))).^(1-σ))).^(1/(1-α*(1-α_m)*(1-σ))),))

            p_o=merge(p_o,(μ_x_c = (1 ./p_o.z_grid_mat) .* ((p_o.Pk/α_m)^α_m) .* ((p_o.w/((1-α)*(1-α_m))).^((1-α)*(1-α_m))).* ((p_o.rtilde_x_c./(α*(1-α_m))).^(α*(1-α_m))),))

            p_o=merge(p_o,(n_x_c = ((1-α)*(1-α_m)/p_o.w)*const_σ*ϕ_x.*(p_o.μ_x_c.^(1-σ)),
            m_x_c = (α_m/p_o.Pk)*const_σ*ϕ_x.*(p_o.μ_x_c.^(1-σ)),

            yd_x_c = const_σ*ϕ_h*(p_o.μ_x_c.^(-σ)),
            yf_x_c = const_σ*(p_o.ξ^σ)*Yf*((1+τ_x)^(-σ))*(τ^(-σ))*(p_o.μ_x_c.^(-σ)),

            pd_x_c = σ/(σ-1).*p_o.μ_x_c,
            pf_x_c = σ/(σ-1)*τ/p_o.ξ.*p_o.μ_x_c))

            #Solution
                p_o.const_x .= p_o.k_const .< real(p_o.k_x_u)
                p_o.k_x[p_o.const_x] .= p_o.k_const[p_o.const_x]
                p_o.n_x[p_o.const_x] .= p_o.n_x_c[p_o.const_x]
                p_o.m_x[p_o.const_x] .= p_o.m_x_c[p_o.const_x]
                p_o.yd_x[p_o.const_x] .= p_o.yd_x_c[p_o.const_x]
                p_o.yf_x[p_o.const_x] .= p_o.yf_x_c[p_o.const_x]
                p_o.pd_x[p_o.const_x] .= p_o.pd_x_c[p_o.const_x]
                p_o.pf_x[p_o.const_x] .= p_o.pf_x_c[p_o.const_x]

        end


    #Profits
        p_o=merge(p_o,(π_x = p_o.pd_x.*p_o.yd_x .+ p_o.ξ*p_o.pf_x.*p_o.yf_x .- ((r+δ)+(1-δ)*cap_gain)*p_o.Pk_lag*p_o.k_x .- p_o.w*p_o.n_x .- p_o.w*F_base .- p_o.Pk*p_o.m_x,))

    ## Non-exporters

        ϕ_nx = ϕ_h

    #Unconstrained
        p_o=merge(p_o,(μ_u = (1 ./p_o.z_grid_mat) * ((p_o.Pk/α_m)^α_m) * ((p_o.w/((1-α)*(1-α_m))).^((1-α)*(1-α_m))) * ((p_o.rtilde_u/(α*(1-α_m))).^(α*(1-α_m))),))

        p_o=merge(p_o,(k_nx_u = (α*(1-α_m)/p_o.rtilde_u)*const_σ*ϕ_nx*(p_o.μ_u.^(1-σ)),
        n_nx_u = ((1-α)*(1-α_m)/p_o.w)*const_σ*ϕ_nx*(p_o.μ_u.^(1-σ)),
        m_nx_u = (α_m/p_o.Pk)*const_σ*ϕ_nx*(p_o.μ_u.^(1-σ)),

        yd_nx_u = const_σ*ϕ_h*(p_o.μ_u.^(-σ))))
        p_o=merge(p_o,(yf_nx_u = zeros(size(p_o.yd_nx_u)),

        pd_nx_u = σ/(σ-1)*p_o.μ_u))
        p_o=merge(p_o,(pf_nx_u = zeros(size(p_o.pd_nx_u)),))

        #Solution
        p_o=merge(p_o,(k_nx = p_o.k_nx_u,
        n_nx = p_o.n_nx_u,
        m_nx = p_o.m_nx_u,
        yd_nx = p_o.yd_nx_u,
        yf_nx = p_o.yf_nx_u,
        pd_nx = p_o.pd_nx_u,
        pf_nx = p_o.pf_nx_u,
        const_nx = BitArray(undef,size(p_o.k_nx_u))))

    #Constrained
        if θ<1+r
            p_o=merge(p_o,(rtilde_nx_c = (((1 ./p_o.k_const) * α * (1-α_m) * const_σ * ϕ_nx) .* ((p_o.μ_u * (p_o.rtilde_u).^(-α*(1-α_m))).^(1-σ))).^(1 ./(1-α*(1-α_m)*(1-σ))),))

            p_o=merge(p_o,(μ_nx_c = (1 ./p_o.z_grid_mat) .* ((p_o.Pk/α_m)^α_m) .* ((p_o.w/((1-α)*(1-α_m))).^((1-α)*(1-α_m))) .* ((p_o.rtilde_nx_c./(α*(1-α_m))).^(α*(1-α_m))),))

            p_o=merge(p_o,(n_nx_c = ((1-α)*(1-α_m)/p_o.w)*const_σ*ϕ_nx*(p_o.μ_nx_c.^(1-σ)),
            m_nx_c = (α_m/p_o.Pk)*const_σ*ϕ_nx*(p_o.μ_nx_c.^(1-σ)),

            yd_nx_c = const_σ*ϕ_h*(p_o.μ_nx_c.^(-σ))))
            p_o=merge(p_o,(yf_nx_c = zeros(size(p_o.yd_nx_c)),

            pd_nx_c = σ/(σ-1)*p_o.μ_nx_c))
            p_o=merge(p_o,(pf_nx_c = zeros(size(p_o.pd_nx_c)),))

            #Solution
            p_o.const_nx .= p_o.k_const .< real(p_o.k_nx_u)
            p_o.k_nx[p_o.const_nx] .= p_o.k_const[p_o.const_nx]
            p_o.n_nx[p_o.const_nx] .= p_o.n_nx_c[p_o.const_nx]
            p_o.m_nx[p_o.const_nx] .= p_o.m_nx_c[p_o.const_nx]
            p_o.yd_nx[p_o.const_nx] .= p_o.yd_nx_c[p_o.const_nx]
            p_o.yf_nx[p_o.const_nx] .= p_o.yf_nx_c[p_o.const_nx]
            p_o.pd_nx[p_o.const_nx] .= p_o.pd_nx_c[p_o.const_nx]
            p_o.pf_nx[p_o.const_nx] .= p_o.pf_nx_c[p_o.const_nx]
        end



    #Profits
        p_o=merge(p_o,(π_nx = p_o.pd_nx.*p_o.yd_nx .- ((r+δ)+(1-δ)*cap_gain)*p_o.Pk_lag*p_o.k_nx .- p_o.w*p_o.n_nx .- p_o.Pk*p_o.m_nx,))

    ## Export decision
        p_o=merge(p_o,(e=BitArray(undef,size(p_o.k_nx_u)),))
        p_o.e .= real(p_o.π_x) .>= real(p_o.π_nx)

      p_o=merge(p_o,(k_const = nothing, n_x_c = nothing, m_x_c = nothing, yd_x_c = nothing, yf_x_c = nothing, pd_x_c = nothing,  pf_x_c = nothing))
      p_o=merge(p_o,(k_x_u = nothing,   n_x_u = nothing, m_x_u = nothing, yd_x_u = nothing, yf_x_u = nothing, pd_x_u = nothing, pf_x_u = nothing))
      p_o=merge(p_o,(n_nx_c = nothing, m_nx_c = nothing, yd_nx_c = nothing, yf_nx_c = nothing, pd_nx_c = nothing, pf_nx_c = nothing))
      p_o=merge(p_o,(k_nx_u = nothing, n_nx_u = nothing, m_nx_u = nothing, yd_nx_u = nothing, yf_nx_u = nothing, pd_nx_u = nothing, pf_nx_u = nothing))

    return p_o

end

########################### Dynamic Problem ########################

function KLS3_dynamicproblem(p_o,p,s,guessV)

@unpack z_P,a_grid,γ,z_grid_size,β,a_grid_size,r=p #See comment at the beginning of KL3_measure

     ### Dynamic problem: value functions and policy functions


     ###Initialize solution objects
        v_new = copy(guessV) #p_o.π_nx
        ap = zeros(Float64,size(v_new))
        ap_ind = zeros(Int64,size(v_new))
        c = similar(ap)

        v_old = copy(v_new)
        v_diff = 1
        exponentg=1-γ

        #Value function iteration algorithm

        profits = p_o.w .+ p_o.tariffsincome .+ p_o.e.*p_o.π_x + (1 .-p_o.e).*p_o.π_nx

        iter = 0
        while v_diff>s.eps
            iter +=1
            v_p=copy(v_old')
            for j = 1:z_grid_size
                @views v_pp = β* z_P[j,:]' *v_p
                Threads.@threads for i = 1:a_grid_size #Threads.@threads is akin to parallelization (in Julia solvers don't have an automatic feature to turn on parallelization, you have to do it manually on the relevant loops [note that you can't do it on all loops, some of them are nor "ordered" in the sense that one iterations builds on the previous one. Those cannot be parallelized])
                    @inbounds begin # @inbounds allows Julia not to check the bounds are correct, thus speeding up the code (however, if bounds are in fact incorrect, this may induce error or throw weird results)
                    c1 =  profits[i,j] .+ a_grid[i].*(1+r) .- a_grid
                    neg_c1_indexes = c1.<=0 #indices of a' for which consumption<0
                    u = (c1.^exponentg)./exponentg .+ neg_c1_indexes.*-1e50

                     v,index_ap = findmax(u .+ v_pp)

                     v_new[i,j] = v
                     ap_ind[i,j]=index_ap[2]
                    end
                end
            end
            Threads.@threads for j =1:z_grid_size
                @inbounds begin
                        @views ap[:,j]= a_grid[ap_ind[:,j]]
                end
            end
            # Accelerator

            # Consumption and utility
             c = profits .+ a_grid'.*ones(1,z_grid_size).*(1+r) .- ap  # The key is ap[:,:], the optimal asset choice

            u = c.^exponentg./exponentg

            for g=1:s.ac_iter
                Threads.@threads for j = 1:z_grid_size
                    @inbounds begin
                    @views v_new[:,j] = u[:,j] .+ (β * z_P[j,:]'*v_new[ap_ind[:,j],:]')'
                    end
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
        p_o=merge(p_o,(v = v_new,
        ap = ap))
        p_o=merge(p_o,(ap_ind=ap_ind,
        c=c,

    ## Store output to be used in simulation

        pd = (1 .-p_o.e).*p_o.pd_nx + p_o.e.*p_o.pd_x,
        yd = (1 .-p_o.e).*p_o.yd_nx + p_o.e.*p_o.yd_x,
        pf = p_o.e.*p_o.pf_x,
        yf = p_o.e.*p_o.yf_x,
        k = (1 .-p_o.e).*p_o.k_nx + p_o.e.*p_o.k_x,
        n = (1 .-p_o.e).*p_o.n_nx + p_o.e.*p_o.n_x,
        m = (1 .-p_o.e).*p_o.m_nx + p_o.e.*p_o.m_x,

        pd_nx = nothing, pd_x = nothing, yd_nx=nothing,  yd_x=nothing, pf_x=nothing, yf_x=nothing, k_nx =nothing, k_x = nothing, n_nx=nothing, n_x=nothing, m_nx=nothing, m_x=nothing,

        π_1 = (1 .-p_o.e).*p_o.π_nx + p_o.e.*p_o.π_x,
        a = p_o.a_grid_mat,

        #Fixed costs
        S_mat = p_o.e*p_o.F,
        F_mat = p_o.e*p_o.F))

    return p_o
end

 ############################## Simulation ##########################
function KLS3_simulate(p_o,p,s,guessM)

@unpack δ,ω_h_k,ω_m_k,ω_m_c,ω_h_c,τ_m_k,τ_m_c,Pm_k,Pm_c,σ,Yf,r=p #See comment at the beginning of KL3_measure
    sim=(measure=KLS3_measure(p_o,p,s,guessM),
    w=p_o.w,
    ξ=p_o.ξ)

    # Exporters and non-exporters
    sim=merge(sim,(mass_x = sum(sim.measure.*p_o.e),))
    sim=merge(sim,(share_x=sim.mass_x, # share of exporters as a fraction of all entrepreneurs
    mass_d = sum(sim.measure.*(1 .-p_o.e)), #Mass of tradable producers who are non-exporters

    # Capital good

    I = δ*sum(sim.measure.*p_o.k),
    M = sum(sim.measure.*p_o.m)))

    yd_k = ((p_o.pd/(p_o.Pk*ω_h_k)).^(-σ)) * (sim.I + sim.M)
    ym_k = ((sim.ξ*Pm_k*(1+τ_m_k)/(p_o.Pk*ω_m_k)).^(-σ)) * (sim.I + sim.M)

    sim=merge(sim,(Yk = (sum(sim.measure.*ω_h_k.*(yd_k.^((σ-1)/σ)) ) + ω_m_k*(ym_k.^((σ-1)/σ)) )^(σ/(σ-1)),
    Pk = (sum(sim.measure.*(ω_h_k^σ).*(p_o.pd.^(1-σ))) + (ω_m_k^σ)*((sim.ξ*Pm_k*(1+τ_m_k)).^(1-σ)) )^(1/(1-σ)),

    # Consumption good

    C = sum(sim.measure.*p_o.c)))

    #yd_c = ((p_o.pd/ω_h_c).^(-σ)) * sim.C
    yd_c = max.(p_o.yd - yd_k,0.00001)
    ym_c = ((sim.ξ*Pm_c*(1+τ_m_c)/ω_m_c).^(-σ)) * sim.C

    sim=merge(sim,(Yc = (sum(sim.measure.*ω_h_c.*(yd_c.^((σ-1)/σ))) + ω_m_c*(ym_c^((σ-1)/σ)))^(σ/(σ-1)),))


    # ϕ_h

    sim=merge(sim,(ϕ_h = (ω_h_c^σ)*sim.Yc + ((ω_h_k*p_o.Pk)^σ)*sim.Yk,

    ### Compute price and quantity indexes

    #Domestic sales in units of final goods
     PdYd =sum(sim.measure[p_o.pd.>0].*p_o.pd[p_o.pd.>0].*p_o.yd[p_o.pd.>0]),

     #Exports
     PxYx = sim.ξ*sum(sim.measure[p_o.pf.>0].*p_o.pf[p_o.pf.>0].*p_o.yf[p_o.pf.>0])))
     sim=merge(sim,(PxYx_USD = sim.PxYx/sim.ξ,

     #Imports
    PmcYmc = sim.ξ*(1+τ_m_c)*Pm_c*ym_c,
    PmkYmk =  sim.ξ*(1+τ_m_k)*Pm_k*ym_k))
    sim=merge(sim,(PmYm = sim.PmcYmc + sim.PmkYmk,

    # Tariffs Income (T)
    tariffsincome=sim.ξ*τ_m_c*Pm_c*ym_c + sim.ξ*τ_m_k*Pm_k*ym_k,

    ### Compute aggregate variables needed to evaluate market clearing conditions

    #Labor and capital
    K = sum(sim.measure.*(p_o.k)))) #recall that already multiplied by oc_choice and sec_choice
    sim=merge(sim,(inv_agg = δ*sim.K,
    N = sum(sim.measure.*(p_o.n))))  # total labor demand for production

    # Fixed costs already corrected by sim.w if needed
    sim=merge(sim,(FC =  sum(sim.measure.*p_o.e.*p_o.F_mat),))

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
    a_supply = sum(sim.measure.*p_o.a_grid_mat),
    a_demand = sim.K))
    #sim=merge(sim,(mc_a = ((1+sim.a_demand)/(1+sim.a_supply)) - 1, #10*log.((1+sim.a_demand)/(1+sim.a_supply)),))
    sim=merge(sim,(mc_a = log.(sim.a_demand/sim.a_supply),

    #Final goods
    #mc_y = ((1+sim.y_demand)/(1+sim.y_supply))-1, #10*log.((1+sim.y_demand)/(1+sim.y_supply)),
    y_supply = sim.Yc,
    y_demand = sim.C))
    sim=merge(sim,(mc_y = log.(sim.y_demand/sim.y_supply),

    k_supply = sim.Yk,
    k_demand = sim.I + sim.M + sim.FC*s.fcost_fgoods))
    sim=merge(sim,(mc_k = log.(sim.k_demand/sim.k_supply),))



    #Beliefs
    #To solve the entrepreneur's problem, they need to have a belief about the aggregate price and quantity indexes in the market
    #In equilibrium, these beliefs need to be consistent with the actual aggregate prices and quantities
    #sim=merge(sim,(mc_y_belief = ((1+sim.ϕ_h)/(1+p_o.ϕ_h))-1,))
    #

    if s.tariffsincome==1
        sim=merge(sim,(mc_Yk_belief = log.(sim.Yk/p_o.Yk),
        mc_tariffs_belief = log.(sim.tariffsincome/p_o.tariffsincome),
        mc_y_belief = ((1+sim.Yc)/(1+p_o.Yc))-1))

    else
        sim=merge(sim,(mc_Yk_belief = 0,
        mc_tariffs_belief = 0,
        mc_y_belief = log.(sim.ϕ_h/p_o.ϕ_h)))
    end

    ## Main statistics

    sim=merge(sim,(Sales = sim.PdYd + sim.PxYx,)) #sim.PxYx is already denominated in units of the domestic final good
    sim=merge(sim,(GDP = sim.Sales - sim.M*p_o.Pk, # Value Added
    gross_output_ratio = (sim.ϕ_h)/(p_o.ξ*Yf),
    gross_output_ratio_2 = (p_o.Yc+p_o.Pk*p_o.Yk)/(p_o.ξ*Yf)))
    sim=merge(sim,(gross_output_ratio_3 =  sim.GDP/(sim.ξ*Yf),

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
    p_o=merge(p_o,(d = (1+r)*(p_o.Pk*p_o.k .- p_o.a),))
    sim=merge(sim,(credit = sum(sim.measure.*max.(p_o.d,0)),)) # only entrepreneurs/firms (for workers p_o.d is negative)
    sim=merge(sim,(credit_gdp = sim.credit/sim.GDP,
    d_agg = sum(sim.measure.*p_o.d))) # for both firms and workers
    sim=merge(sim,(NFA_GDP = -sim.d_agg/sim.GDP,

    credit_sales = sim.credit/(sim.PdYd + sim.PxYx),

    k_wagebill = p_o.Pk*sim.K/(sim.w*sim.n_supply)))

    sales_x =  sim.ξ*p_o.pf.*p_o.yf #p_opf is denominated in foreign currency, so we adjust it
    sales_d = p_o.pd.*p_o.yd
    sales = sales_d+sales_x
    x_share = sales_x./sales

    sim=merge(sim,(x_share_av =  sum(filter(!isnan,sim.measure.*p_o.e.* x_share))./ sum(filter(!isnan,sim.measure.*p_o.e)),))
    #sim.x_share_av[isnan.(sim.x_share_av)].=1  #Don't know the role of this, but it throws error in Julia given sim.x_share_av is a scalar


    SalesExporters = sum(sim.measure.*sales.*p_o.e)
    sim=merge(sim,(x_share_wav =  Base.sum(filter(!isnan,sim.measure.*p_o.e.*(sales./SalesExporters).*x_share))/ sum(sim.measure.*p_o.e.*(sales./SalesExporters)),))
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

    sim=merge(sim,(sales_avg_nx = sum(sim.measure.*(1 .-p_o.e).*sales) / sum(sim.measure.*(1 .-p_o.e)),
    labor_avg_nx = sum(sim.measure.*(1 .-p_o.e).*(p_o.n )) / sum(sim.measure.*(1 .-p_o.e))))
    sim=merge(sim,(sales_d_avg_nx = sim.sales_avg_nx,

    sales_avg_x = sum(sim.measure.*p_o.e.*sales)/ sum(sim.measure.*p_o.e),
    labor_avg_x = ((sim.FC/sim.w)*(1-s.fcost_fgoods) + sum(sim.measure.*p_o.e.*(p_o.n))) / sum(sim.measure.*p_o.e),

    # These include fixed costs
    labor_tot = sim.N,  # total labor demand for production
    labor_tot_x = (sim.FC/sim.w)*(1-s.fcost_fgoods)  + sum(sim.measure.*p_o.e.*(p_o.n)),
    labor_tot_nx =  sum(sim.measure.*(1 .-p_o.e).*(p_o.n)),
    labor_tot_w = sum(sim.measure)))

    sim=merge(sim,(labor_x_share = sim.labor_tot_x/sim.labor_tot,
    labor_nx_share = sim.labor_tot_nx/sim.labor_tot,

    sales_d_avg_x = sum(sim.measure.*p_o.e.*sales_d)/ sum(sim.measure.*p_o.e),

    xpremium_sales = sim.sales_avg_x/sim.sales_avg_nx,
    xpremium_labor = sim.labor_avg_x/sim.labor_avg_nx))

    sim=merge(sim,(xpremium_sales_d = sim.sales_d_avg_x/sim.sales_d_avg_nx,

    ln_sales_sd_mean = sim.ln_sales_sd /  sim.ln_sales_mean,
    ln_sales_d_sd_mean = sim.ln_sales_d_sd /  sim.ln_sales_d_mean,
    sales_sd_mean = sim.sales_sd /  sales_mean,
    sales_d_sd_mean = sim.sales_d_sd /  sales_d_mean))

    labor_mean= sum(sim.measure.*(p_o.n+ (p_o.e.*p_o.F_mat./sim.w )*(1-s.fcost_fgoods)))
    labor_sd=sqrt.(sum(sim.measure.*(p_o.n+ (p_o.e.*p_o.F_mat./sim.w)*(1-s.fcost_fgoods) .- labor_mean).^2))
    sim=merge(sim,(labor_sd_mean = labor_sd /  labor_mean,


   # Average productivity
    avg_productivity1 = ((sum(sim.measure.*(p_o.z_grid_mat).^(σ-1)))/sum(sim.measure)).^(1/(σ - 1)),
    avg_productivity2 =  sum(sim.measure.*sales.* p_o.z_grid_mat) / sum(sim.measure .* sales),
    avg_productivity3 =  sum(sim.measure.*sales_d.* p_o.z_grid_mat) / sum(sim.measure .* sales_d),



 ## Solution-related statistics

# All
 a_min_share = sum(sim.measure[1,:]),
 a_max_share = sum(sim.measure[end,:])))

 return sim, p_o

end
############################# GE_par (2 versions, first for using NLsolve [in-place function], second standard)#############################

####### For NLsolve #######
#This version of KLS3_GE_par is in-place, that is, its output (market clearing conditions) is a modified value of an input (F, in this case). This is done because NLsolve needs this format to work. It is customary to add a ! at the end of in-place functions.
function KLS3_GE_par!(F,x,p,s)

    @unpack ω_h_c, ω_h_k,ω_m_c, ω_m_k, σ, Pm_c, Pm_k, τ_m_c, τ_m_k, F_base,z_π,a_grid_size,z_grid_size,r,tariffsincome=p #See comment at the beginning of KL3_measure

    #### Aggregate prices and quantities ####
    #Guessed prices

    p_o=(w = exp(x[1]),
    ϕ_h = exp(x[2]),
    ξ = exp(x[3]),
    Pk = exp(x[4]),
    Pk_lag=exp(x[4]))


    ### Solution

    # Tariff income (initialized to 0)
    if s.tariffsincome == 1
        p_o=(w = exp(x[1]),
        ϕ_h = exp(x[2]),
        ξ = exp(x[3]),
        Pk = exp(x[4]),
        Pk_lag=exp(x[4]),
        Yk = exp(x[5]),
        Yc = exp(x[2]))
        p_o=merge(p_o,(ϕ_h = (ω_h_c^σ)*p_o.Yc + ((ω_h_k*p_o.Pk)^σ)*p_o.Yk,))

        #Yc =  (p_o.ϕ_h - ( ω_m_k*p_o.Pk)^σ*p_o.Yk)/( ω_m_c^σ)
        Yc = p_o.Yc
        ym_c = Yc*(p_o.ξ*Pm_c*(1+τ_m_c)/ω_m_c)^(-σ)
        ym_k = p_o.Yk*(p_o.Pk^σ)*(p_o.ξ*Pm_k*(1+τ_m_k)/ω_m_k)^(-σ)

        if s.tariffsincome_PE==0
            p_o=merge(p_o,(tariffsincome = τ_m_c*p_o.ξ*Pm_c*ym_c + τ_m_k*p_o.ξ*Pm_k*ym_k,))
        else
            p_o=merge(p_o,(tariffsincome=p.tariffsincome,))
        end

    end

    #Display
     if s.display==1
         print("\nGuesses: ϕ_h=$(p_o.ϕ_h) , w= $(p_o.w) , r= $(r) , ξ= $(p_o.ξ) , Pk= $(p_o.Pk)")
     end


     #Fixed costs
     if s.fcost_fgoods==0 # If in units of labor
         p_o=merge(p_o,(F = p_o.w*F_base,))
      else
         p_o=merge(p_o,(F = F_base,))
     end

    #Static problems for tradable and nontradable sectors
    p_o_static=KLS3_staticproblem(p_o,p)
    p_o=merge(p_o,p_o_static)


    #Dynamic problem and simulation (No sunk costs)
    global guessV=copy(p_o.π_nx)

    p_o_dynamic=KLS3_dynamicproblem(p_o,p,s,guessV)
    p_o=merge(p_o,p_o_dynamic)

    global guessM=vcat(z_π',zeros(a_grid_size-1,z_grid_size))

    sim,p_o1 = KLS3_simulate(p_o,p,s,guessM)
    p_o=merge(p_o,p_o1)

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
if s.display==1 #\n is for advancing the display one line downward (i.e. to display each printed string a different line)

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

    return F #This is important for NLsolve to work.

end

####### Standard #######

function KLS3_GE_par(x,p,s)
#See comments of the in-place version (KLS3_GE_par!)
@unpack ω_h_c, ω_h_k,ω_m_c, ω_m_k, σ, Pm_c, Pm_k, τ_m_c, τ_m_k, F_base,z_π,a_grid_size,z_grid_size,r=p

#### Aggregate prices and quantities ####
#Guessed prices

p_o=(w = exp(x[1]),
ϕ_h = exp(x[2]),
ξ = exp(x[3]),
Pk = exp(x[4]),
Pk_lag=exp(x[4]))


### Solution

# Tariff income (initialized to 0)
if s.tariffsincome == 1
    p_o=(w = exp(x[1]),
    ϕ_h = exp(x[2]),
    ξ = exp(x[3]),
    Pk = exp(x[4]),
    Pk_lag=exp(x[4]),
    Yk = exp(x[5]),
    Yc = exp(x[2]))
    p_o=merge(p_o,(ϕ_h = (ω_h_c^σ)*p_o.Yc + ((ω_h_k*p_o.Pk)^σ)*p_o.Yk,))

    #Yc =  (p_o.ϕ_h - ( ω_m_k*p_o.Pk)^σ*p_o.Yk)/( ω_m_c^σ)
    Yc = p_o.Yc
    ym_c = Yc*(p_o.ξ*Pm_c*(1+τ_m_c)/ω_m_c)^(-σ)
    ym_k = p_o.Yk*(p_o.Pk^σ)*(p_o.ξ*Pm_k*(1+τ_m_k)/ω_m_k)^(-σ)

    if s.tariffsincome_PE==0
        p_o=merge(p_o,(tariffsincome = τ_m_c*p_o.ξ*Pm_c*ym_c + τ_m_k*p_o.ξ*Pm_k*ym_k,))
    else
        p_o=merge(p_o,(tariffsincome=p.tariffsincome,))
    end

end

#Display
 if s.display==1
     print("\nGuesses: ϕ_h=$(p_o.ϕ_h) , w= $(p_o.w) , r= $(r) , ξ= $(p_o.ξ) , Pk= $(p_o.Pk)")
 end


 #Fixed costs
 if s.fcost_fgoods==0 # If in units of labor
     p_o=merge(p_o,(F = p_o.w*F_base,))
  else
     p_o=merge(p_o,(F = p_o.F_base,))
 end

 #Static problems for tradable and nontradable sectors
     p_o_static=KLS3_staticproblem(p_o,p)
     p_o=merge(p_o,p_o_static)

guessV=copy(p_o.π_nx)

#Dynamic problem and simulation (No sunk costs)
p_o_dynamic=KLS3_dynamicproblem(p_o,p,s,guessV)
p_o=merge(p_o,p_o_dynamic)

guessM=vcat(z_π',zeros(a_grid_size-1,z_grid_size))

sim,p_o1 = KLS3_simulate(p_o,p,s,guessM)
p_o=merge(p_o,p_o1)


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

 return mcc, p_o, sim

end
