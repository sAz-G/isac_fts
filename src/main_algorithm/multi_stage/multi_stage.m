
%%
clear;
clc;
close all;

addpath(genpath("..\..\..\src"));


%% Simulation parameter

run("call_hyperParam.m")


eta = params.sim.eta; 
mu  = params.sim.mu; 
M = 5;
N_tot = 80;
N_stg = 25;
K_stg = floor(N_stg/mu);
K_tot = floor(N_tot/mu);

% Energy parameters
E_total = 40e3; % [J]
E_m = E_total;
E_min = 5e3;

%% Simulation Setup
% Basestation
S_s = [100; 100];
% Communication user
S_c = [1.3e3; 1.2e3];
% Sensing target
S_t = [200; 1.3e3];

% estimate the target randomly 
S_target_est = S_t + [100; -100];

% Middle point between communication user und target user
S_mid = (S_c + S_target_est)/2;

% Inital trajectory
[S_init, V_init] = init_trajectory(S_s, N_stg, S_mid, V_str, T_f);

plot_map(S_init, S_s, S_t, S_target_est, S_c);




%% Optimization
w_star  = params.sim.w_star;
iter    = params.sim.iter;
S_opt_mat = nan(2,N_stg,M);
m = 1;

D_meas = nan(K_stg,M);

while E_min < E_m

% Middle point between communication user und target user
S_mid = (S_c + S_target_est)/2;

% Inital trajectory
[S_init, V_init] = init_trajectory(S_s, N_stg, S_mid, V_str, T_f);

% calc initial values of parametes
delta_square_init = sqrt(1 + norms(V_init, 2, 1).^4/(4*v_0^4)) - norms(V_init, 2 ,1).^2/(2 * v_0^2);
xi_init = delta_square_init;

% run the mth stag
[S_opt_m,E_m_used] = single_stage(E_m, N_stg, delta_square_init,K_stg, S_c, S_init,S_target_est,S_s,V_init,params);

for k = 1:K_stg
D_meas(k,m) = sense_target(S_t, S_opt_m(:,k*mu));
end

% get the new estimated target
[x_t,y_t] = grid_vectors(1500,1500,1000,1000);

[x_t_idx,y_t_idx]  = get_min(D_meas(:,m),x_t,y_t,S_opt_m(1, mu:mu:N_stg),S_opt_m(2, mu:mu:N_stg), H,params);

S_target_est = [x_t(x_t_idx), y_t(y_t_idx)]; 
% calculate the energy 
E_m = E_m - E_m_used; 

% store calculated trajectory 
S_opt_mat(:,1:size(S_opt_m,2), m) = S_opt_m;

% increase the iteration variable
m = m+1;

end

% last stage here
plot_map(S_opt_m, S_s, S_t, S_target_est, S_c);

