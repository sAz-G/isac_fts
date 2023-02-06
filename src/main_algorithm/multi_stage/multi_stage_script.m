
%%
clear;
clc;
close all;

addpath(genpath("..\..\..\src"));


%% Simulation parameter

run("call_hyperParam.m")

params.sim.iter = 20;
params.sim.eta = 0.5;

mu  = params.sim.mu; 
M = 5;
N_tot = nan;
N_stg = params.sim.N_stg;
K_stg = floor(N_stg/mu);
K_tot = nan;

% Energy parameters
E_total = 30e3; % [J]
E_m = E_total;
E_min = calc_real_energy(K_stg, params.sim.V_max*ones(1,N_stg), params);

%% Simulation Setup
% Basestation
S_b = [100; 100];
S_s = S_b;
% Communication user
S_c = [1.3e3; 1.2e3];
% Sensing target
S_t = [200; 1.3e3];

% estimate the target randomly 
S_target_est = S_t + [100; -100];

%% Optimization
m = 1;

D_meas           = nan(K_stg,M);
S_opt_mat        = nan(2,N_stg,M);
S_init_mat       = nan(2,N_stg,M);
S_target_est_mat = nan(2,M);
E_min_vec        = nan(M,1);
E_m_vec          = nan(M,1);
E_used_vec       = nan(M,1);
V_m_mat          = nan(2,N_stg,M);
xi_m_vec         = nan(N_stg,1,M);
delta_m_vec      = nan(N_stg,1,M);
R_opt_vecs       = nan(params.sim.iter, 1, M);
CRB_opt_vecs     = nan(params.sim.iter, 1, M);

while E_min < E_m

% Middle point between communication user und target user
S_mid = (S_c + S_target_est)/2;
% S_mid = S_target_est*params.sim.eta + S_c*(1-params.sim.eta);

% Inital trajectory
[S_init, V_init] = init_trajectory(S_s, S_mid, N_stg, params);

plot_map(S_init, S_b, S_t, S_target_est, S_c, params);

% calc initial values of parametes
delta_square_init = sqrt(1 + norms(V_init, 2, 1).^4/(4*params.energy.v_0^4)) - norms(V_init, 2 ,1).^2/(2*params.energy.v_0^2);

% run the mth stage
debug_mode = true;
if debug_mode == false
    [S_opt_m,E_m_used, V_m, xi_m, delta_m,CRB_vec_m,R_vec_m] = single_stage(E_m, N_stg, delta_square_init,K_stg, S_c, S_init,S_target_est,S_s,V_init,params);
else
    [S_opt_m, V_m, xi_m, delta_m,CRB_vec_m,R_vec_m] = single_stage_debug(E_m, N_stg, delta_square_init, K_stg, S_c, S_init, S_target_est, S_s, V_init, m,params);
end

D_meas(:,m) = sense_target(S_t, S_opt_m(:,mu:mu:end), params);

% store calculated trajectory 
S_opt_mat(:,1:size(S_opt_m,2), m) = S_opt_m;

V_froms_S = calc_velocity_from_trajectory(S_opt_m, S_s, params);
E_m_used = calc_real_energy(K_stg, V_froms_S, params);

% get the new estimated target
[x_t,y_t] = grid_vectors(1500,1500,1000,1000);

[x_t_idx,y_t_idx]  = get_min(D_meas(1:K_stg*m),x_t,y_t,reshape(S_opt_mat(1, mu:mu:end,1:m),1,m*K_stg), ...
                                                       reshape(S_opt_mat(2, mu:mu:end,1:m),1,m*K_stg),params);
                                                
S_target_est = [x_t(x_t_idx); y_t(y_t_idx)]; 

S_init_mat(:,1:size(S_opt_m,2), m) = S_init;
% set new starting point
S_s = S_opt_m(:,end);

% calculate the energy 
E_m = E_m - E_m_used; 
E_min = calc_back_energy(S_opt_m(:,end), S_b, params);

%store variables
E_used_vec(m)           = E_m_used;
E_m_vec(m)              = E_m;
E_min_vec(m)            = E_min;
S_target_est_mat(:,m)   = S_target_est;
V_m_mat(:,:,m)          = V_m;
delta_m_vec(:,:,m)      = delta_m;
xi_m_vec(:,:,m)         = xi_m;
CRB_opt_vecs(:,:,m)     = CRB_vec_m;
R_opt_vecs(:,:,m)       = R_vec_m;

% increase the iteration variable
m = m+1;

end

M = m-1;
N_tot = size(S_opt_mat,2)*size(S_opt_mat,3);
K_tot = floor(N_tot/mu);

% last stage here
plot_map(S_opt_mat, S_b, S_t, S_target_est, S_c,params);

