
%%
clear;
clc;
close all;

addpath(genpath("..\..\..\src"));


%% Simulation parameter
run("call_hyperParam.m")

mu  = params.sim.mu;      % measurment step. Measure sensing target every mu steps
M = 7;                      
N_stg = params.sim.N_stg; % amount of flight points in every stage 
K_stg = floor(N_stg/mu);  % amount of hovering points in every stage 

%% Simulation Setup
% Basestation
s_b = setup.base_station_pos;        % base station position 
s_s = s_b;                           % start position of the initial trajectory (not necessary of the real trajectory)

% Communication user
s_c = setup.comm_user_pos;           % position of the communication user 
% Sensing target
s_t = setup.sense_target_pos;        % real position of the sensing target

% estimate the target randomly 
s_target_est = setup.est_sense_target; % target estimation at the beginning 

% Energy parameters
E_total = setup.total_energy;   % total amount of energy in [J]  
E_m = E_total;                  % amount of energy at the mth stage, start with E_total 

% the farthest point the quad can fly to from its current position 
s_max_end = s_s + [1/sqrt(2) * params.sim.T_f * params.sim.V_max * N_stg; 1/sqrt(2) * params.sim.T_f * params.sim.V_max * N_stg];
E_min = calc_back_energy(s_max_end, s_b, params); % amount of energy needed to fly back to the bas statio 

%% variables to save 
m = 1;  % first stage 
D_meas           = nan(K_stg,M);                   % measurement matrix
S_opt_mat        = nan(2,N_stg,M);                 % matrix of optimal traj. for each stage
S_init_mat       = nan(2,N_stg,M);                 % matrix of initial trajectories
S_target_est_mat = nan(2,M);                       % matrix of target estimations of each stage 
E_min_vec        = nan(M,1);                       % vector that stores the minimum energy needed at each stage 
E_m_vec          = nan(M,1);                       % store the energy availale at each stage 
E_used_vec       = nan(M,1);                       % store the energy used at each stage 
V_m_mat          = nan(2,N_stg,M);                 % store the optimal velocity at each stage
xi_m_vec         = nan(N_stg,1,M);                 % store the optimal xi at each stage 
delta_m_vec      = nan(N_stg,1,M);                 % store the optimal delta at each stage
R_opt_vecs       = nan(params.sim.iter, M);        % store optimal data rate solutions 
CRB_opt_vecs     = nan(params.sim.iter, M);        % store the optimalcrb soutions 

%% Optimization
epsilon = params.sim.eta;
if params.sim.eta == 1
    epsilon = 0.99;
end

S_total_m   = [];
S_hover_total     = [];

while E_min < E_m
 
        
% end point between communication user und target user for the initial traj
s_end = s_target_est*epsilon + s_c*(1-epsilon);

% Save the whole trajectory. Use a [2, (m-1)*N_stg] structure 
%S_total_m = [S_total_m, S_opt_m]; %reshape(S_opt_mat(:,:,1:m-1), [2, (m-1)*N_stg]); % total trajectory at the mth stage

% Inital trajectory
[S_init, V_init] = init_trajectory(s_s, s_end, N_stg, params);
plot_map(S_init, s_b, s_t, s_target_est, s_c, params);

% get hover points
S_total_m = horzcat(S_total_m,S_init);
[~,S_hover_total] = get_S_hover(params, 1, m, S_total_m);

% calc initial values of parametes
delta_square_init = sqrt(1 + norms(V_init, 2, 1).^4/(4*params.energy.v_0^4))...
                            - norms(V_init, 2 ,1).^2/(2*params.energy.v_0^2);

% run the mth stage
debug_mode = true;

if debug_mode == false
    [S_opt_m, E_m_used, V_m, xi_m, delta_m,CRB_vec_m,R_vec_m] = optimize_m(E_m, N_stg, delta_square_init,K_stg, s_c, S_init,s_target_est,s_s,V_init,params);
else
    [S_opt_m, V_m, xi_m, delta_m,CRB_vec_m,R_vec_m] = optimize_m_debug(E_m, s_c, S_hover_total, S_total_m,delta_square_init, s_target_est, s_s, V_init, params);
end

% sense target at each hover point 
D_meas(:,m) = sense_target(s_t, S_opt_m(:,mu:mu:end), params);

% store calculated trajectory 
S_opt_mat(:,1:size(S_opt_m,2), m) = S_opt_m;

% calculate the used energy
E_m_used = calc_real_energy(K_stg, S_opt_m, s_s, params);

% get the new estimated target
[x_t,y_t] = grid_vectors(1500,1500,1000,1000); % create vectors for grid search

% get the estimated target matrix index
[x_t_idx,y_t_idx]  = get_min(D_meas(1:K_stg*m),x_t,y_t,...
                             reshape(S_opt_mat(1, mu:mu:end,1:m),1,m*K_stg), ...
                             reshape(S_opt_mat(2, mu:mu:end,1:m),1,m*K_stg),...
                             params);
                                                
s_target_est = [x_t(x_t_idx); y_t(y_t_idx)]; 

S_init_mat(:,1:size(S_opt_m,2), m) = S_init;
% set new current point
s_s = S_opt_m(:,end);

% Save the whole trajectory. Use a [2, (m-1)*N_stg] structure 
S_total_m(:,end-N_stg+1:end) = S_opt_m;

% calculate the energy 
E_m = E_m - E_m_used; 
E_min = calc_back_energy(S_opt_m(:,end), s_b, params);

%store variables
E_used_vec(m)           = E_m_used;
E_m_vec(m)              = E_m;
E_min_vec(m)            = E_min;
S_target_est_mat(:,m)   = s_target_est;
V_m_mat(:,:,m)          = V_m;
delta_m_vec(:,:,m)      = delta_m;
xi_m_vec(:,:,m)         = xi_m;
CRB_opt_vecs(:,m)       = CRB_vec_m;
R_opt_vecs(:,m)         = R_vec_m;

% increase the iteration variable
m = m+1;

end

M = m-1;                                     % amount of stages 
N_tot = size(S_opt_mat,2)*size(S_opt_mat,3); % total amount of points 
K_tot = floor(N_tot/mu);                     % total amount of hover points 

% last stage here
plot_map(S_opt_mat, s_b, s_t, s_target_est, s_c,params);  % plot map 

