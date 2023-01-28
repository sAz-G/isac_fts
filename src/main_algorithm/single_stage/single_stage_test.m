%%
clear;
clc;
close all;

addpath(genpath(".\crb_functions"));
addpath(genpath(".\helper_functions"));
addpath(genpath(".\rate_functions"));

%% Simulation parameter
run("call_hyperParam.m")
eta = .5; 
mu  = 5; 

%% Simulation Setup
% Basestation
S_s = [100; 100];
% Communication user
S_c = [1.3e3; 1.2e3];
% Sensing target
S_t = [200; 1.3e3];

N_tot = 25;

%% Multi-stage approach for UWV trajectory design - First stage

% Number of setpoints for the drone
K_tot = floor(N_tot/mu); % Number of hover points

% Energy parameters
E_total = 40e3; % [J]
E_m = E_total;

% target estimation via grid search
S_target_est = S_t ;%+ [100; -100];

% Middle point between communication user und target user
S_mid = (S_c + S_target_est)/2;

% Inital trajectory
[S_init, V_init] = init_trajectory(S_s, N_tot, S_mid, V_str, T_f);
plot_map(S_init, S_s, S_t, S_target_est, S_c);

delta_square_last = sqrt(1 + norms(V_init, 2, 1).^4/(4*v_0^4)) - norms(V_init, 2 ,1).^2/(2 * v_0^2);


%% Optimization
w_star = 0.2;
iter = 5;


S_opt = single_stage(H,eta, E_m, N_tot, delta_square_last,w_star, K_tot,iter, mu,S_c, S_init,S_target_est,S_s,V_init);
plot_map(S_opt, S_s, S_t, S_target_est, S_c);
