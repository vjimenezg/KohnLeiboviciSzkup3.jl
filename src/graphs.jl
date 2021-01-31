####################### Test Graphs #######################

T=s.N;

#Real GDP
plot(1:T,[0 (sim_fun.GDP_Laspeyres[2:T]./sim_0.GDP .-1)']',title= "GDP Laspeyres baseline weighted 1 tau")


#Capital
plot(1:T,[0 (sim_fun.K[2:T]./sim_0.K.-1)']', title="K baseline weighted 1 tau")

#Consumption
plot(1:T,[0 (sim_fun.C[2:T]./sim_0.C.-1)']',title="C baseline weighted 1 tau")

#Real Exports
plot(1:T,[0 (sim_fun.X_Laspeyres[2:T]./sim_0.PxYx.-1)']',title="X Laspeyres baseline weighted 1 tau")

#Real Domestic sales
plot(1:T,[0 (sim_fun.D_Laspeyres[2:T]./sim_0.PdYd.-1)']',title="D Laspeyres baseline weighted 1 tau")

#Real Exchange Rate
plot(1:T,[0 (sim_fun.ξ[2:T]./sim_0.ξ.-1)']',title="RER baseline weighted 1 tau")

#Real Wage
plot(1:T,[0 (sim_fun.w[2:T]./sim_0.w.-1)']',title=
"Wage baseline weighted 1 tau")

#Price of Capital
plot(1:T,[0 (sim_fun.Pk[2:T]./sim_0.Pk.-1)']',title="Pk baseline weighted 1 tau")

#Phi_h
plot(1:T,[0 (sim_fun.ϕ_h[2:T]./sim_0.ϕ_h.-1)']',title="Phih baseline weighted 1 tau")
