## Welfare analysis

q=(v0_transition = rt[2].v,
q.v0_transition_t3 = rt[3].v,
q.v0_ss = rt[1].v,
q.v0_ss_final = rt[s.N].v,

#Individual-specific gains
ind_welfare = (q.v0_transition./q.v0_ss).^(1/(1-m.γ)) .- 1,
ind_welfare_t3 = (q.v0_transition_t3./q.v0_ss).^(1/(1-m.γ)) .- 1,
measure_0 = sim_0.measure,
measure_end = sim_end.measure,
ind_welfare_avg = sum(q.measure_0.*q.ind_welfare),

#Gains without transition
ind_welfare_notrans = (q.v0_ss_final./q.v0_ss).^(1/(1-m.γ)) .- 1,
ind_welfare_avg_notrans = sum(q.measure_0.*q.ind_welfare_notrans),

#Average gains
agg_welfare_avg = (sum(q.measure_0.*q.v0_transition)/sum(q.measure_0.*q.v0_ss)).^(1/(1-m.γ)) .- 1,

#Gains without transition
agg_welfare_avg_notrans = (sum(q.measure_end.*q.v0_ss_final)/sum(q.measure_0.*q.v0_ss)).^(1/(1-m.γ)) .- 1,
agg_welfare_avg_notrans_fixedmeasure = (sum(q.measure_0.*q.v0_ss_final)/sum(q.measure_0.*q.v0_ss)).^(1/(1-m.γ)) .- 1,
agg_welfare_avg_notrans_fixedvalue = (sum(q.measure_end.*q.v0_ss)/sum(q.measure_0.*q.v0_ss)).^(1/(1-m.γ)) .- 1)


### STATISTICS

# All
q=merge(q,(share_winners = sum(sim_0.measure[q.ind_welfare.>0]),
share_losers = sum(sim_0.measure[q.ind_welfare.<0]),
welfare_avg_winners = sum(sim_0.measure[q.ind_welfare.>=0].*q.ind_welfare[q.ind_welfare.>=0])/ sum(sim_0.measure[q.ind_welfare.>=0]),
welfare_avg_losers = sum(sim_0.measure[q.ind_welfare.<0].*q.ind_welfare[q.ind_welfare.<0])/ sum(sim_0.measure[q.ind_welfare.<0])))

# Exporters vs Non-exporters
q=merge(q,(welfare_avg_exp = sum(sim_0.measure[rt[1].e.==1].*q.ind_welfare[rt[1].e.==1])/sum(sim_0.measure[rt[1].e.==1]),
welfare_avg_nonexp = sum(sim_0.measure[rt[1].e.==0].*q.ind_welfare[rt[1].e.==0])/sum(sim_0.measure[rt[1].e.==0]),

share_losers_exp = sum(sim_0.measure[q.ind_welfare.<0 .& rt[1].e.==1])/ sum(sim_0.measure[rt[1].e.==1]),
share_losers_nonexp = sum(sim_0.measure[q.ind_welfare.<0 .& rt[1].e.==0])/  sum(sim_0.measure[rt[1].e.==0]),
share_winners_exp = sum(sim_0.measure[q.ind_welfare.>=0 .& rt[1].e.==1]) / sum(sim_0.measure[rt[1].e.==1]),
share_winners_nonexp = sum(sim_0.measure[q.ind_welfare.>=0 .& rt[1].e.==0]) / sum(sim_0.measure[rt[1].e.==0]),

welfare_avg_losers_exp=  sum(sim_0.measure[rt[1].e.==1].*q.ind_welfare[rt[1].e.==1].*(q.ind_welfare[rt[1].e.==1].<0))/sum(sim_0.measure[rt[1].e.==1].*(q.ind_welfare[rt[1].e.==1].<0)),
welfare_avg_winners_exp=  sum(sim_0.measure[rt[1].e.==1].*q.ind_welfare[rt[1].e.==1].*(q.ind_welfare[rt[1].e.==1].>=0))/sum(sim_0.measure[rt[1].e.==1].*(q.ind_welfare[rt[1].e.==1].>=0)),
welfare_avg_losers_nonexp=  sum(sim_0.measure[rt[1].e.==0].*q.ind_welfare[rt[1].e.==0].*(q.ind_welfare[rt[1].e.==0].<0))/sum(sim_0.measure[rt[1].e.==0].*(q.ind_welfare[rt[1].e.==0].<0)),
welfare_avg_winners_nonexp=  sum(sim_0.measure[rt[1].e.==0].*q.ind_welfare[rt[1].e.==0].*(q.ind_welfare[rt[1].e.==0].>=0))/sum(sim_0.measure[rt[1].e.==0].*(q.ind_welfare[rt[1].e.==0].>=0))))


# Rich vs Poor
q=merge(q,(measure_a_cum = cumsum(sum(sim_0.measure,dims=2),dims=1),
aP = sum(q.measure_a_cum.<0.99,dims=1) .+1,
# aP = s.a_grid_size./2 .+1; # threshold for asset

welfare_avg_poor = sum(sim_0.measure[1:q.aP-1,:].*q.ind_welfare[1:q.aP-1,:])/sum(sim_0.measure[1:q.aP-1,:]),
welfare_avg_rich = sum(sim_0.measure[q.aP:end,:].*q.ind_welfare[q.aP:end,:])/sum(sim_0.measure[q.aP:end,:]),

share_losers_poor = sum(sim_0.measure[1:q.aP-1,:].*(q.ind_welfare[1:q.aP-1,:].<0))/sum(sim_0.measure[1:q.aP-1,:]),
share_winners_poor = sum(sim_0.measure[1:q.aP-1,:].*(q.ind_welfare[1:q.aP-1,:].>=0))/sum(sim_0.measure[1:q.aP-1,:]),
welfare_avg_losers_poor =  sum(sim_0.measure[1:q.aP-1,:].*q.ind_welfare[1:q.aP-1,:].*(q.ind_welfare[1:q.aP-1,:].<0))/sum(sim_0.measure[1:q.aP-1,:].*(q.ind_welfare[1:q.aP-1,:].<0)),
welfare_avg_winners_poor =  sum(sim_0.measure[1:q.aP-1,:].*q.ind_welfare[1:q.aP-1,:].*(q.ind_welfare[1:q.aP-1,:].>=0))/sum(sim_0.measure[1:q.aP-1,:].*(q.ind_welfare[1:q.aP-1,:].>=0)),

share_losers_rich = sum(sim_0.measure[q.aP:end,:].*(q.ind_welfare[q.aP:end,:].<0))/sum(sim_0.measure[q.aP:end,:]),
share_winners_rich = sum(sim_0.measure[q.aP:end,:].*(q.ind_welfare[q.aP:end,:].>=0))/sum(sim_0.measure[q.aP:end,:]),
welfare_avg_losers_rich = sum(sim_0.measure[q.aP:end,:].*q.ind_welfare[q.aP:end,:].*(q.ind_welfare[q.aP:end,:].<0))/sum(sim_0.measure[q.aP:end,:].*(q.ind_welfare[q.aP:end,:].<0)),
welfare_avg_winners_rich = sum(sim_0.measure[q.aP:end,:].*q.ind_welfare[q.aP:end,:].*(q.ind_welfare[q.aP:end,:].>=0))/sum(sim_0.measure[q.aP:end,:].*(q.ind_welfare[q.aP:end,:].>=0))))

# Productive vs Non-productive
q=merge(q,(measure_z_cum = cumsum(sum(sim_0.measure,dims=1),dims=2),
zP = sum(q.measure_z_cum.<0.5,dims=1).+1,
# zP = s.z_grid_size./2 .+1; # threshold for productive

welfare_avg_NP = sum(sim_0.measure[:,1:q.zP-1].*q.ind_welfare[:,1:q.zP-1])/sum(sim_0.measure[:,1:q.zP-1]),
welfare_avg_P = sum(sim_0.measure[:,q.zP:end].*q.ind_welfare[:,q.zP:end])/sum(sim_0.measure[:,q.zP:end]),
share_losers_NP =  sum(sim_0.measure[:,1:q.zP-1].*(q.ind_welfare[:,1:q.zP-1].<0))/sum(sim_0.measure[:,1:q.zP-1]),
share_winners_NP =  sum(sim_0.measure[:,1:q.zP-1].*(q.ind_welfare[:,1:q.zP-1].>=0))/sum(sim_0.measure[:,1:q.zP-1]),
welfare_avg_winners_NP  = sum(sim_0.measure[:,1:q.zP-1].*q.ind_welfare[:,1:q.zP-1].*(q.ind_welfare[:,1:q.zP-1].>=0))/sum(sim_0.measure[:,1:q.zP-1].*(q.ind_welfare[:,1:q.zP-1].>=0)),
welfare_avg_losers_NP  = sum(sim_0.measure[:,1:q.zP-1].*q.ind_welfare[:,1:q.zP-1].*(q.ind_welfare[:,1:q.zP-1].<0))/sum(sim_0.measure[:,1:q.zP-1].*(q.ind_welfare[:,1:q.zP-1].<0)),
share_winners_P=  sum(sim_0.measure[:,q.zP:end].*(q.ind_welfare[:,q.zP:end].>=0))/sum(sim_0.measure[:,q.zP:end]),
share_losers_P =  sum(sim_0.measure[:,q.zP:end].*(q.ind_welfare[:,q.zP:end].<0))/sum(sim_0.measure[:,q.zP:end]),
welfare_avg_losers_P = sum(sim_0.measure[:,q.zP:end].*q.ind_welfare[:,q.zP:end].*(q.ind_welfare[:,q.zP:end].<0))/sum(sim_0.measure[:,q.zP:end].*(q.ind_welfare[:,q.zP:end].<0)),
welfare_avg_winners_P = sum(sim_0.measure[:,q.zP:end].*q.ind_welfare[:,q.zP:end].*(q.ind_welfare[:,q.zP:end].>=0))/sum(sim_0.measure[:,q.zP:end].*(q.ind_welfare[:,q.zP:end].>=0))))


###  Exporter vs Nonexporters and Rich vs Poor
q=merge(q,(measure_a_cum = cumsum(sum(sim_0.measure,dims=2),dims=1),
aP = sum(q.measure_a_cum .<0.99,dims=1) .+1,
# aP = s.a_grid_size/2+1, # threshold for asset
aVP = 1,

welfare_avg_expR =  sum((rt[1].e[q.aP:end,:].==1).*sim_0.measure[q.aP:end,:].*q.ind_welfare[q.aP:end,:])/sum((rt[1].e[q.aP:end,:].==1).*sim_0.measure[q.aP:end,:]),
welfare_avg_nonexpR =  sum((rt[1].e[q.aP:end,:].==0).*sim_0.measure[q.aP:end,:].*q.ind_welfare[q.aP:end,:])/sum((rt[1].e[q.aP:end,:].==0).*sim_0.measure[q.aP:end,:]),
welfare_avg_expP =  sum((rt[1].e[q.aVP+1:q.aP-1,:].==1).*sim_0.measure[q.aVP+1:q.aP-1,:].*q.ind_welfare[q.aVP+1:q.aP-1,:])/sum((rt[1].e[q.aVP+1:q.aP-1,:].==1).*sim_0.measure[q.aVP+1:q.aP-1,:]),
welfare_avg_nonexpP =  sum((rt[1].e[q.aVP+1:q.aP-1,:].==0).*sim_0.measure[q.aVP+1:q.aP-1,:].*q.ind_welfare[q.aVP+1:q.aP-1,:])/sum((rt[1].e[q.aVP+1:q.aP-1,:].==0).*sim_0.measure[q.aVP+1:q.aP-1,:]),
welfare_avg_expVP = sum((rt[1].e[1:q.aVP ,:].==1).*sim_0.measure[1:q.aVP ,:].*q.ind_welfare[1:q.aVP ,:])/sum((rt[1].e[1:q.aVP ,:].==1).*sim_0.measure[1:q.aVP ,:]),
welfare_avg_nonexpVP = sum((rt[1].e[1:q.aVP ,:].==0).*sim_0.measure[1:q.aVP ,:].*q.ind_welfare[1:q.aVP ,:])/sum((rt[1].e[1:q.aVP ,:].==0).*sim_0.measure[1:q.aVP ,:]),


share_expR = sum((rt[1].e[q.aP:end,:].==1).*sim_0.measure[q.aP:end,:]),
share_expP = sum((rt[1].e[q.aVP+1:q.aP-1,:].==1).*sim_0.measure[q.aVP+1:q.aP-1,:]),
share_nonexpR = sum((rt[1].e[q.aP:end,:].==0).*sim_0.measure[q.aP:end,:]),
share_nonexpP = sum((rt[1].e[q.aVP+1:q.aP-1,:].==0).*sim_0.measure[q.aVP+1:q.aP-1,:]),
share_expVP = sum((rt[1].e[1:q.aVP,:].==1).*sim_0.measure[1:q.aVP,:]),
share_nonexpVP = sum((rt[1].e[1:q.aVP,:].==0).*sim_0.measure[1:q.aVP,:]),


share_losers_expR = sum((q.ind_welfare[q.aP:end,:].<0) .*  (rt[1].e[q.aP:end,:].==1).*sim_0.measure[q.aP:end,:]),
share_winners_expR = sum((q.ind_welfare[q.aP:end,:].>=0) .*  (rt[1].e[q.aP:end,:].==1).*sim_0.measure[q.aP:end,:]),
share_losers_nonexpR = sum( (q.ind_welfare[q.aP:end,:].<0) .*  (rt[1].e[q.aP:end,:].==0).*sim_0.measure[q.aP:end,:]),
share_winners_nonexpR = sum( (q.ind_welfare[q.aP:end,:].>=0) .*  (rt[1].e[q.aP:end,:].==0).*sim_0.measure[q.aP:end,:]),

share_losers_expP = sum( (q.ind_welfare[q.aVP+1:q.aP-1,:].<0) .*  (rt[1].e[q.aVP+1:q.aP-1,:].==1).*sim_0.measure[q.aVP+1:q.aP-1,:]),
share_winners_expP = sum( (q.ind_welfare[q.aVP+1:q.aP-1,:].>=0) .*  (rt[1].e[q.aVP+1:q.aP-1,:].==1).*sim_0.measure[q.aVP+1:q.aP-1,:]),
share_losers_nonexpP = sum( (q.ind_welfare[q.aVP+1:q.aP-1,:].<0) .*  (rt[1].e[q.aVP+1:q.aP-1,:].==0).*sim_0.measure[q.aVP+1:q.aP-1,:]),
share_winners_nonexpP = sum((q.ind_welfare[q.aVP+1:q.aP-1,:].>=0) .*  (rt[1].e[q.aVP+1:q.aP-1,:].==0).*sim_0.measure[q.aVP+1:q.aP-1,:]),

share_losers_expVP = sum((q.ind_welfare[1:q.aVP,:].<0) .*  (rt[1].e[1:q.aVP,:].==1).*sim_0.measure[1:q.aVP,:]),
share_winners_expVP = sum((q.ind_welfare[1:q.aVP,:].>=0) .*  (rt[1].e[1:q.aVP,:].==1).*sim_0.measure[1:q.aVP,:]),
share_losers_nonexpVP = sum((q.ind_welfare[1:q.aVP,:].<0) .*  (rt[1].e[1:q.aVP,:].==0).*sim_0.measure[1:q.aVP,:]),
share_winners_nonexpVP = sum((q.ind_welfare[1:q.aVP,:].>=0) .*  (rt[1].e[1:q.aVP,:].==0).*sim_0.measure[1:q.aVP,:]),



welfare_avg_losers_expR = sum(sim_0.measure[q.aP:end,:] .* q.ind_welfare[q.aP:end,:] .* (q.ind_welfare[q.aP:end,:].<0) .*(rt[1].e[q.aP:end,:].==1))sum(sim_0.measure[q.aP:end,:].*(q.ind_welfare[q.aP:end,:].<0) .*  (rt[1].e[q.aP:end,:].==1)),
welfare_avg_winners_expR = sum(sim_0.measure[q.aP:end,:].*q.ind_welfare[q.aP:end,:].*(q.ind_welfare[q.aP:end,:].>=0) .* (rt[1].e[q.aP:end,:].==1))/sum(sim_0.measure[q.aP:end,:].*(q.ind_welfare[q.aP:end,:].>=0) .* (rt[1].e[q.aP:end,:].==1)),
welfare_avg_losers_nonexpR =sum(sim_0.measure[q.aP:end,:] .* q.ind_welfare[q.aP:end,:] .* (q.ind_welfare[q.aP:end,:].<0) .*  (rt[1].e[q.aP:end,:].==0))/sum(sim_0.measure[q.aP:end,:].*(q.ind_welfare[q.aP:end,:].<0) .*  (rt[1].e[q.aP:end,:].==0)),
welfare_avg_winners_nonexpR = sum(sim_0.measure[q.aP:end,:].*q.ind_welfare[q.aP:end,:].*(q.ind_welfare[q.aP:end,:].>=0) .* (rt[1].e[q.aP:end,:].==0))/sum(sim_0.measure[q.aP:end,:].*(q.ind_welfare[q.aP:end,:].>=0).*(rt[1].e[q.aP:end,:].==0)),


welfare_avg_losers_expP = sum(sim_0.measure[q.aVP+1:q.aP-1,:] .*q.ind_welfare[q.aVP+1:q.aP-1,:].* (q.ind_welfare[q.aVP+1:q.aP-1,:].<0) .*  (rt[1].e[q.aVP+1:q.aP-1,:].==1))/sum(sim_0.measure[q.aVP+1:q.aP-1,:].*(q.ind_welfare[q.aVP+1:q.aP-1,:].<0) .*  (rt[1].e[q.aVP+1:q.aP-1,:].==1) ),
welfare_avg_winners_expP = sum(sim_0.measure[q.aVP+1:q.aP-1,:].*q.ind_welfare[q.aVP+1:q.aP-1,:].*(q.ind_welfare[q.aVP+1:q.aP-1,:].>=0) .* (rt[1].e[q.aVP+1:q.aP-1,:].==1))/sum(sim_0.measure[q.aVP+1:q.aP-1,:].*(q.ind_welfare[q.aVP+1:q.aP-1,:].>=0) .* (rt[1].e[q.aVP+1:q.aP-1,:].==1) ),
welfare_avg_losers_nonexpP = sum(sim_0.measure[q.aVP+1:q.aP-1,:] .* q.ind_welfare[q.aVP+1:q.aP-1,:] .* (q.ind_welfare[q.aVP+1:q.aP-1,:].<0) .*  (rt[1].e[q.aVP+1:q.aP-1,:].==0))/sum(sim_0.measure[q.aVP+1:q.aP-1,:].*(q.ind_welfare[q.aVP+1:q.aP-1,:].<0) .*  (rt[1].e[q.aVP+1:q.aP-1,:].==0) ),
welfare_avg_winners_nonexpP = sum(sim_0.measure[q.aVP+1:q.aP-1,:].*q.ind_welfare[q.aVP+1:q.aP-1,:].*(q.ind_welfare[q.aVP+1:q.aP-1,:].>=0) .* (rt[1].e[q.aVP+1:q.aP-1,:].==0))/sum(sim_0.measure[q.aVP+1:q.aP-1,:].*(q.ind_welfare[q.aVP+1:q.aP-1,:].>=0) .* (rt[1].e[q.aVP+1:q.aP-1,:].==0) ),


welfare_avg_losers_expVP = sum(sim_0.measure[1:q.aVP,:] .* q.ind_welfare[1:q.aVP,:] .* (q.ind_welfare[1:q.aVP,:].<0) .*  (rt[1].e[1:q.aVP,:].==1))/sum(sim_0.measure[1:q.aVP,:].*(q.ind_welfare[1:q.aVP,:].<0) .*  (rt[1].e[1:q.aVP,:].==1)),
q.welfare_avg_winners_expVP = sum(sim_0.measure[1:q.aVP,:].*q.ind_welfare[1:q.aVP,:].*(q.ind_welfare[1:q.aVP,:].>=0) .* (rt[1].e[1:q.aVP,:].==1))/sum(sim_0.measure[1:q.aVP,:].*(q.ind_welfare[1:q.aVP,:].>=0) .* (rt[1].e[1:q.aVP,:].==1)),
q.welfare_avg_losers_nonexpVP =sum(sim_0.measure[1:q.aVP,:] .* q.ind_welfare[1:q.aVP,:] .* (q.ind_welfare[1:q.aVP,:].<0) .*  (rt[1].e[1:q.aVP,:].==0))/sum(sim_0.measure[1:q.aVP,:].*(q.ind_welfare[1:q.aVP,:].<0) .*  (rt[1].e[1:q.aVP,:].==0)),
q.welfare_avg_winners_nonexpVP = sum(sim_0.measure[1:q.aVP,:].*q.ind_welfare[1:q.aVP,:].*(q.ind_welfare[1:q.aVP,:].>=0) .* (rt[1].e[1:q.aVP,:].==0))/sum(sim_0.measure[1:q.aVP,:].*(q.ind_welfare[1:q.aVP,:].>=0) .* (rt[1].e[1:q.aVP,:].==0))))



###  Productive vs Non-productive and Rich vs Poor
q=merge(q,(zP_low = sum(q.measure_z_cum .<0.33,dims=1) .+1,
zP_med = sum(q.measure_z_cum.<0.66,dims=1) .+1,

# shares
share_lowR = sum(sim_0.measure[q.aP:end,1:q.zP_low.-1]),
share_lowP = sum(sim_0.measure[q.aVP.+1:q.aP-1,1:q.zP_low.-1]),
share_lowVP = sum(sim_0.measure[1:q.aVP,1:q.zP_low.-1]),

share_medR = sum(sim_0.measure[q.aP:end,q.zP_low:q.zP_med.-1]),
share_medP = sum(sim_0.measure[q.aVP.+1:q.aP.-1,q.zP_low:q.zP_med.-1]),
share_medVP = sum(sim_0.measure[1:q.aVP,q.zP_low:q.zP_med]),

share_highR = sum(sim_0.measure[q.aP:end,q.zP_med:end]),
share_highP = sum(sim_0.measure[q.aVP.+1:q.aP.-1,q.zP_med:end]),
share_highVP = sum(sim_0.measure[1:q.aVP,q.zP_med:end])))


# welfare
q=merge(q,(welfare_avg_lowR =  sum(sim_0.measure[q.aP:end,1:q.zP_low.-1].*q.ind_welfare[q.aP:end,1:q.zP_low.-1])/sum(sim_0.measure[q.aP:end,1:q.zP_low.-1]),
welfare_avg_lowP =  sum(sim_0.measure[q.aVP.+1:q.aP.-1,1:q.zP_low.-1].*q.ind_welfare[q.aVP.+1:q.aP.-1,1:q.zP_low.-1])/sum(sim_0.measure[q.aVP.+1:q.aP.-1,1:q.zP_low.-1]),
welfare_avg_lowVP =  sum(sim_0.measure[1:q.aVP ,1:q.zP_low.-1].*q.ind_welfare[1:q.aVP ,1:q.zP_low.-1])/sum(sim_0.measure[1:q.aVP ,1:q.zP_low.-1]),

welfare_avg_medR =  sum(sim_0.measure[q.aP:end,q.zP_low:q.zP_med.-1].*q.ind_welfare[q.aP:end,q.zP_low:q.zP_med.-1])/sum(sim_0.measure[q.aP:end,q.zP_low:q.zP_med.-1]),
welfare_avg_medP =  sum(sim_0.measure[q.aVP.+1:q.aP.-1,q.zP_low:q.zP_med.-1].*q.ind_welfare[q.aVP.+1:q.aP.-1,q.zP_low:q.zP_med.-1])/sum(sim_0.measure[q.aVP.+1:q.aP.-1,q.zP_low:q.zP_med.-1]),
welfare_avg_medVP =  sum(sim_0.measure[1:q.aVP ,q.zP_low:q.zP_med.-1].*q.ind_welfare[1:q.aVP ,q.zP_low:q.zP_med.-1])/sum(sim_0.measure[1:q.aVP ,q.zP_low:q.zP_med.-1]),

welfare_avg_highR =  sum(sim_0.measure[q.aP:end,q.zP_med:end].*q.ind_welfare[q.aP:end,q.zP_med:end])/sum(sim_0.measure[q.aP:end,q.zP_med:end]),
welfare_avg_highP =  sum(sim_0.measure[q.aVP.+1:q.aP.-1,q.zP_med:end].*q.ind_welfare[q.aVP.+1:q.aP.-1,q.zP_med:end])/sum(sim_0.measure[q.aVP.+1:q.aP.-1,q.zP_med:end]),
welfare_avg_highVP =  sum(sim_0.measure[1:q.aVP ,q.zP_med:end].*q.ind_welfare[1:q.aVP ,q.zP_med:end])/sum(sim_0.measure[1:q.aVP ,q.zP_med:end])))

# share of losers
q=merge(q,(share_losers_lowR = sum((q.ind_welfare[q.aP:end,1:q.zP_low.-1].<0) .* sim_0.measure[q.aP:end,1:q.zP_low.-1]),
share_losers_lowP = sum((q.ind_welfare[q.aVP.+1:q.aP.-1,1:q.zP_low.-1].<0) .* sim_0.measure[q.aVP.+1:q.aP.-1,1:q.zP_low.-1]),
share_losers_lowVP = sum((q.ind_welfare[1:q.aVP,1:q.zP_low.-1].<0) .*  sim_0.measure[1:q.aVP,1:q.zP_low-.1]),

share_losers_medR = sum((q.ind_welfare[q.aP:end,q.zP_low:q.zP_med.-1].<0) .*  sim_0.measure[q.aP:end,q.zP_low:q.zP_med.-1]),
share_losers_medP = sum((q.ind_welfare[q.aVP.+1:q.aP.-1,q.zP_low:q.zP_med.-1].<0) .*sim_0.measure[q.aVP.+1:q.aP.-1,q.zP_low:q.zP_med.-1]),
share_losers_medVP = sum((q.ind_welfare[1:q.aVP,q.zP_low:q.zP_med.-1].<0) .*sim_0.measure[1:q.aVP,q.zP_low:q.zP_med.-1]),

share_losers_highR = sum((q.ind_welfare[q.aP:end,q.zP_med:end].<0) .*  sim_0.measure[q.aP:end,q.zP_med:end]),
share_losers_highP = sum((q.ind_welfare[q.aVP+1:q.aP-1,q.zP_med:end].<0) .* sim_0.measure[q.aVP.+1:q.aP.-1,q.zP_med:end]),
share_losers_highVP = sum((q.ind_welfare[1:q.aVP,q.zP_med:end].<0) .* sim_0.measure[1:q.aVP,q.zP_med:end])))


# share of winners
q=merge(q,(share_winners_lowR = sum((q.ind_welfare[q.aP:end,1:q.zP_low.-1].>=0) .* sim_0.measure[q.aP:end,1:q.zP_low.-1]),
share_winners_lowP = sum((q.ind_welfare[q.aVP.+1:q.aP.-1,1:q.zP_low.-1].>=0) .* sim_0.measure[q.aVP.+1:q.aP.-1,1:q.zP_low.-1]),
share_winners_lowVP = sum((q.ind_welfare[1:q.aVP,1:q.zP_low.-1].>=0) .*  sim_0.measure[1:q.aVP,1:q.zP_low.-1]),

share_winners_medR = sum((q.ind_welfare[q.aP:end,q.zP_low:q.zP_med.-1].>=0) .*  sim_0.measure[q.aP:end,q.zP_low:q.zP_med.-1]),
share_winners_medP = sum((q.ind_welfare[q.aVP.+1:q.aP.-1,q.zP_low:q.zP_med.-1].>=0) .*sim_0.measure[q.aVP.+1:q.aP.-1,q.zP_low:q.zP_med.-1]),
share_winners_medVP = sum((q.ind_welfare[1:q.aVP,q.zP_low:q.zP_med.-1].>=0) .*sim_0.measure[1:q.aVP,q.zP_low:q.zP_med.-1]),

share_winners_highR = sum((q.ind_welfare[q.aP:end,q.zP_med:end].>=0) .*  sim_0.measure[q.aP:end,q.zP_med:end]),
share_winners_highP = sum((q.ind_welfare[q.aVP.+1:q.aP.-1,q.zP_med:end].>=0) .* sim_0.measure[q.aVP.+1:q.aP.-1,q.zP_med:end]),
share_winners_highVP = sum((q.ind_welfare[1:q.aVP,q.zP_med:end].>=0) .* sim_0.measure[1:q.aVP,q.zP_med:end])))

# avg welfare of losers
q=merge(q,(welfare_avg_losers_lowR = sum(sim_0.measure[q.aP:end,1:q.zP_low-1] .* q.ind_welfare[q.aP:end,1:q.zP_low-1] .* (q.ind_welfare[q.aP:end,1:q.zP_low-1].<0))/sum(sim_0.measure[q.aP:end,1:q.zP_low-1].*(q.ind_welfare[q.aP:end,1:q.zP_low-1].<0)),
welfare_avg_losers_lowP =sum(sim_0.measure[q.aVP+1:q.aP-1,1:q.zP_low-1] .* q.ind_welfare[q.aVP+1:q.aP-1,1:q.zP_low-1] .* (q.ind_welfare[q.aVP+1:q.aP-1,1:q.zP_low-1].<0))/sum(sim_0.measure[q.aVP+1:q.aP-1,1:q.zP_low-1].*(q.ind_welfare[q.aVP+1:q.aP-1,1:q.zP_low-1].<0)),
welfare_avg_losers_lowVP = sum(sim_0.measure[1:q.aVP,1:q.zP_low-1] .* q.ind_welfare[1:q.aVP,1:q.zP_low-1] .* (q.ind_welfare[1:q.aVP,1:q.zP_low-1].<0))/sum(sim_0.measure[1:q.aVP,1:q.zP_low-1].*(q.ind_welfare[1:q.aVP,1:q.zP_low-1].<0)),


welfare_avg_losers_medR = sum(sim_0.measure[q.aP:end,q.zP_low:q.zP_med.-1] .* q.ind_welfare[q.aP:end,q.zP_low:q.zP_med.-1] .* (q.ind_welfare[q.aP:end,q.zP_low:q.zP_med.-1].<0))/sum(sim_0.measure[q.aP:end,q.zP_low:q.zP_med.-1].*(q.ind_welfare[q.aP:end,q.zP_low:q.zP_med.-1].<0)),
welfare_avg_losers_medP = sum(sim_0.measure[q.aVP.+1:q.aP.-1,q.zP_low:q.zP_med.-1] .* q.ind_welfare[q.aVP.+1:q.aP.-1,q.zP_low:q.zP_med.-1] .* (q.ind_welfare[q.aVP.+1:q.aP.-1,q.zP_low:q.zP_med.-1].<0))/sum(sim_0.measure[q.aVP.+1:q.aP.-1,q.zP_low:q.zP_med.-1].*(q.ind_welfare[q.aVP.+1:q.aP.-1,q.zP_low:q.zP_med.-1].<0)),
welfare_avg_losers_medVP = sum(sim_0.measure[1:q.aVP,q.zP_low:q.zP_med.-1] .* q.ind_welfare[1:q.aVP,q.zP_low:q.zP_med.-1] .* (q.ind_welfare[1:q.aVP,q.zP_low:q.zP_med.-1].<0))/sum(sim_0.measure[1:q.aVP,q.zP_low:q.zP_med.-1].*(q.ind_welfare[1:q.aVP,q.zP_low:q.zP_med.-1].<0)),


welfare_avg_losers_highR = sum(sim_0.measure[q.aP:end,q.zP_med:end] .* q.ind_welfare[q.aP:end,q.zP_med:end] .* (q.ind_welfare[q.aP:end,q.zP_med:end].<0))/sum(sim_0.measure[q.aP:end,q.zP_med:end].*(q.ind_welfare[q.aP:end,q.zP_med:end].<0)),
welfare_avg_losers_highP = sum(sim_0.measure[q.aVP.+1:q.aP.-1,q.zP_med:end] .* q.ind_welfare[q.aVP.+1:q.aP.-1,q.zP_med:end] .* (q.ind_welfare[q.aVP.+1:q.aP.-1,q.zP_med:end].<0))/sum(sim_0.measure[q.aVP.+1:q.aP.-1,q.zP_med:end].*(q.ind_welfare[q.aVP.+1:q.aP.-1,q.zP_med:end].<0)),
welfare_avg_losers_highVP = sum(sim_0.measure[1:q.aVP,q.zP_med:end] .* q.ind_welfare[1:q.aVP,q.zP_med:end] .* (q.ind_welfare[1:q.aVP,q.zP_med:end].<0))/sum(sim_0.measure[1:q.aVP,q.zP_med:end].*(q.ind_welfare[1:q.aVP,q.zP_med:end].<0))))

# avg welfare of winners
q=merge(q,(welfare_avg_winners_lowR = sum(sim_0.measure[q.aP:end,1:q.zP_low.-1] .* q.ind_welfare[q.aP:end,1:q.zP_low.-1] .* (q.ind_welfare[q.aP:end,1:q.zP_low.-1].>=0))/sum(sim_0.measure[q.aP:end,1:q.zP_low.-1].*(q.ind_welfare[q.aP:end,1:q.zP_low.-1].>=0)),
welfare_avg_winners_lowP = sum(sim_0.measure[q.aVP.+1:q.aP.-1,1:q.zP_low.-1] .* q.ind_welfare[q.aVP.+1:q.aP.-1,1:q.zP_low.-1] .* (q.ind_welfare[q.aVP.+1:q.aP.-1,1:q.zP_low.-1].>=0))/sum(sim_0.measure[q.aVP.+1:q.aP.-1,1:q.zP_low.-1].*(q.ind_welfare[q.aVP.+1:q.aP.-1,1:q.zP_low.-1].>=0)),
welfare_avg_winners_lowVP = sum(sim_0.measure[1:q.aVP,1:q.zP_low.-1] .* q.ind_welfare[1:q.aVP,1:q.zP_low.-1] .* (q.ind_welfare[1:q.aVP,1:q.zP_low.-1].>=0))/sum(sim_0.measure[1:q.aVP,1:q.zP_low.-1].*(q.ind_welfare[1:q.aVP,1:q.zP_low.-1].>=0)),


welfare_avg_winners_medR = sum(sim_0.measure[q.aP:end,q.zP_low:q.zP_med.-1] .* q.ind_welfare[q.aP:end,q.zP_low:q.zP_med.-1] .* (q.ind_welfare[q.aP:end,q.zP_low:q.zP_med.-1].>=0))/sum(sim_0.measure[q.aP:end,q.zP_low:q.zP_med.-1].*(q.ind_welfare[q.aP:end,q.zP_low:q.zP_med.-1].>=0)),
welfare_avg_winners_medP = sum(sim_0.measure[q.aVP.+1:q.aP.-1,q.zP_low:q.zP_med.-1] .* q.ind_welfare[q.aVP.+1:q.aP.-1,q.zP_low:q.zP_med.-1] .* (q.ind_welfare[q.aVP.+1:q.aP.-1,q.zP_low:q.zP_med.-1].>=0))/sum(sim_0.measure[q.aVP.+1:q.aP.-1,q.zP_low:q.zP_med.-1].*(q.ind_welfare[q.aVP.+1:q.aP.-1,q.zP_low:q.zP_med.-1].>=0)),
welfare_avg_winners_medVP = sum(sim_0.measure[1:q.aVP,q.zP_low:q.zP_med.-1] .* q.ind_welfare[1:q.aVP,q.zP_low:q.zP_med.-1] .* (q.ind_welfare[1:q.aVP,q.zP_low:q.zP_med.-1].>=0))/sum(sim_0.measure[1:q.aVP,q.zP_low:q.zP_med.-1].*(q.ind_welfare[1:q.aVP,q.zP_low:q.zP_med.-1].>=0)),


welfare_avg_winners_highR = sum(sim_0.measure[q.aP:end,q.zP_med:end] .* q.ind_welfare[q.aP:end,q.zP_med:end] .* (q.ind_welfare[q.aP:end,q.zP_med:end].>=0))/sum(sim_0.measure[q.aP:end,q.zP_med:end].*(q.ind_welfare[q.aP:end,q.zP_med:end].>=0)),
welfare_avg_winners_highP = sum(sim_0.measure[q.aVP.+1:q.aP.-1,q.zP_med:end] .* q.ind_welfare[q.aVP.+1:q.aP.-1,q.zP_med:end] .* (q.ind_welfare[q.aVP.+1:q.aP.-1,q.zP_med:end].>=0))/sum(sim_0.measure[q.aVP.+1:q.aP.-1,q.zP_med:end].*(q.ind_welfare[q.aVP.+1:q.aP.-1,q.zP_med:end].>=0)),
welfare_avg_winners_highVP = sum(sim_0.measure[1:q.aVP,q.zP_med:end] .* q.ind_welfare[1:q.aVP,q.zP_med:end] .* (q.ind_welfare[1:q.aVP,q.zP_med:end].>=0))/sum(sim_0.measure[1:q.aVP,q.zP_med:end].*(q.ind_welfare[1:q.aVP,q.zP_med:end].>=0))))


q=merge(q,(table_exp=zeros(6,6),table_prod=zeros(6,6)))

### TABLES EXPORT STATUS & NET WORTH
q.table_exp[1,:] = [q.share_nonexpVP q.share_nonexpP q.share_nonexpR q.share_expVP q.share_expP q.share_expR]
q.table_exp[2,:] = [q.welfare_avg_nonexpVP q.welfare_avg_nonexpP q.welfare_avg_nonexpR q.welfare_avg_expVP q.welfare_avg_expP q.welfare_avg_expR]
q.table_exp[3,:] = [q.share_winners_nonexpVP q.share_winners_nonexpP q.share_winners_nonexpR q.share_winners_expVP q.share_winners_expP q.share_winners_expR]
q.table_exp[4,:] = [q.welfare_avg_winners_nonexpVP q.welfare_avg_winners_nonexpP q.welfare_avg_winners_nonexpR q.welfare_avg_winners_expVP q.welfare_avg_winners_expP q.welfare_avg_winners_expR]
q.table_exp[5,:] = [q.share_losers_nonexpVP q.share_losers_nonexpP q.share_losers_nonexpR q.share_losers_expVP q.share_losers_expP q.share_losers_expR]
q.table_exp[6,:] = [q.welfare_avg_losers_nonexpVP q.welfare_avg_losers_nonexpP q.welfare_avg_losers_nonexpR q.welfare_avg_losers_expVP q.welfare_avg_losers_expP q.welfare_avg_losers_expR]



### TABLES PRODUCTIVITY & NET WORTH

q.table_prod[1,:] = [q.share_lowVP q.share_lowP q.share_lowR q.share_medVP q.share_medP q.share_medR q.share_highVP q.share_highP q.share_highR]
q.table_prod[2,:] = [q.welfare_avg_lowVP q.welfare_avg_lowP q.welfare_avg_lowR q.welfare_avg_medVP q.welfare_avg_medP q.welfare_avg_medR q.welfare_avg_highVP q.welfare_avg_highP q.welfare_avg_highR]
q.table_prod[3,:] = [q.share_winners_lowVP q.share_winners_lowP q.share_winners_lowR q.share_winners_medVP q.share_winners_medP q.share_winners_medR q.share_winners_highVP q.share_winners_highP q.share_winners_highR]
q.table_prod[4,:] = [q.welfare_avg_winners_lowVP q.welfare_avg_winners_lowP q.welfare_avg_winners_lowR q.welfare_avg_winners_medVP q.welfare_avg_winners_medP q.welfare_avg_winners_medR q.welfare_avg_winners_highVP q.welfare_avg_winners_highP q.welfare_avg_winners_highR]
q.table_prod[5,:] = [q.share_losers_lowVP q.share_losers_lowP q.share_losers_lowR q.share_losers_medVP q.share_losers_medP q.share_losers_medR q.share_losers_highVP q.share_losers_highP q.share_losers_highR]
q.table_prod[6,:] = [q.welfare_avg_losers_lowVP q.welfare_avg_losers_lowP q.welfare_avg_losers_lowR q.welfare_avg_losers_medVP q.welfare_avg_losers_medP q.welfare_avg_losers_medR q.welfare_avg_losers_highVP q.welfare_avg_losers_highP q.welfare_avg_losers_highR]


### Bins
sav =0
zbin = 5
abin = 5

q=merge(q,(measure_bins_5=zeros(abin,zbin),measure_bins_5_end=zeros(abin,zbin),avg_welfare_bins_5=zeros(abin,zbin),winners_bins_5=zeros(abin,zbin),losers_bins_5=zeros(abin,zbin)))

for j = 1:zbin
    for i = 1:abin
        pL = (j-1)/zbin*s.z_grid_size+1
        pH = j/zbin*s.z_grid_size
        aL = (i-1)/abin*s.a_grid_size+1
        aH = i/abin*s.a_grid_size

        q.measure_bins_5[i,j] = sum(sim_0.measure[aL:aH,pL:pH]);
        q.measure_bins_5_end[i,j] = sum(sim_end.measure[aL:aH,pL:pH]),
        q.avg_welfare_bins_5[i,j] = sum(sim_0.measure[aL:aH,pL:pH] .* q.ind_welfare[aL:aH,pL:pH])./sum(sim_0.measure[aL:aH,pL:pH]),
        q.winners_bins_5[i,j] = sum(sim_0.measure[aL:aH,pL:pH] .* (q.ind_welfare[aL:aH,pL:pH].>=0))/sum(sim_0.measure[aL:aH,pL:pH]),
        q.losers_bins_5[i,j] = sum(sim_0.measure[aL:aH,pL:pH] .* (q.ind_welfare[aL:aH,pL:pH].<0))/sum(sim_0.measure[aL:aH,pL:pH])
    end
end


#FROM HERE UNTIL "3D bar (10 bins)" PENDING


figure(1)
b1 = bar3(q.measure_bins_5); colorbar
for k = 1:length(b1)
    zdata = b1(k).ZData;
    b1(k).CData = zdata;
    b1(k).FaceColor = 'interp';
end
ylabel('assets')
xlabel('productivity')
# view(32, 22)
if sav == 1
        saveas(gcf,'Figures/measure_FF_tau_c_5','epsc')
        saveas(gcf,'Figures/measure_FF_tau_c_5','pdf')
end

figure(2)
b2 = bar3(q.avg_welfare_bins_5); colorbar
for k = 1:length(b2)
    zdata = b2(k).ZData;
    b2(k).CData = zdata;
    b2(k).FaceColor = 'interp';
end
xlabel('productivity')
ylabel('assets')
# view(32, 22)
if sav == 1
        saveas(gcf,'Figures/welfare_FF_tau_c_5','epsc')
        saveas(gcf,'Figures/welfare_FF_tau_c_5','pdf')
end

figure(3)
b3 = bar3(q.losers_bins_5); colorbar
for k = 1:length(b3)
    zdata = b3(k).ZData;
    b3(k).CData = zdata;
    b3(k).FaceColor = 'interp';
end
xlabel('productivity')
ylabel('assets')
set(gca,'Zlim',[0,1])
caxis([0 1])
# view(32, 22)
if sav == 1
        saveas(gcf,'Figures/welfare_FF_tau_c_5','epsc')
        saveas(gcf,'Figures/welfare_FF_tau_c_5','pdf')
end


figure(4)
imagesc(q.winners_bins_5); colorbar
xlabel('productivity')
ylabel('assets')
# view(32, 22)
if sav == 1
        saveas(gcf,'Figures/winners_FF_tau_c_5','epsc')
        saveas(gcf,'Figures/winners_FF_tau_c_5','pdf')
end


figure(5)
imagesc(q.ind_welfare); colorbar
if sav == 1
        saveas(gcf,'Figures/welfare_image','epsc')
        saveas(gcf,'Figures/welfare_image','pdf')
end

### 3D bar (10 bins)
zbin = 10;
abin = 10;


q=merge(q,(measure_bins_10=zeros(abin,zbin),avg_welfare_bins_10=zeros(abin,zbin),winners_bins_10=zeros(abin,zbin),losers_bins_10=zeros(abin,zbin)))


for j = 1:zbin
    for i = 1:abin
        pL = (j-1)/zbin*s.z_grid_size+1;
        pH = j/zbin*s.z_grid_size;
        aL = (i-1)/abin*s.a_grid_size+1;
        aH = i/abin*s.a_grid_size;

        q.measure_bins_10[i,j] = sum(sim_0.measure[aL:aH,pL:pH])
        q.avg_welfare_bins_10[i,j] = sum(sim_0.measure[aL:aH,pL:pH] .* q.ind_welfare[aL:aH,pL:pH])./sum(sim_0.measure[aL:aH,pL:pH])
        q.winners_bins_10[i,j] = sum(sim_0.measure[aL:aH,pL:pH] .* (q.ind_welfare[aL:aH,pL:pH].>=0))/ sum(sim_0.measure[aL:aH,pL:pH])
        q.losers_bins_10[i,j] = sum(sim_0.measure[aL:aH,pL:pH] .* (q.ind_welfare[aL:aH,pL:pH].<0))/sum(sim_0.measure[aL:aH,pL:pH])
    end
end

#FROM HERE UNTIL THE END PENDING
figure(5)
b1 = bar3(q.measure_bins_10); colorbar
for k = 1:length(b1)
    zdata = b1(k).ZData;
    b1(k).CData = zdata;
    b1(k).FaceColor = 'interp';
end
ylabel('assets')
xlabel('productivity')
# view(32, 22)
if sav == 1
        saveas(gcf,'Figures/measure_FF_tau_c_10','epsc')
        saveas(gcf,'Figures/measure_FF_tau_c_10','pdf')
end

figure(6)
b2 = bar3(q.avg_welfare_bins_10); colorbar
for k = 1:length(b2)
    zdata = b2(k).ZData;
    b2(k).CData = zdata;
    b2(k).FaceColor = 'interp';
end
xlabel('productivity')
ylabel('assets')
# view(32, 22)
if sav == 1
        saveas(gcf,'Figures/welfare_FF_tau_c_10','epsc')
        saveas(gcf,'Figures/welfare_FF_tau_c_10','pdf')
end

figure(7)
imagesc(q.winners_bins_10); colorbar
# xlabel('productivity')
# ylabel('assets')
# view(32, 22)
if sav == 1
        saveas(gcf,'Figures/winners_FF_tau_c_10','epsc')
        saveas(gcf,'Figures/winners_FF_tau_c_10','pdf')
end
