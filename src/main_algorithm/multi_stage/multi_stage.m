function results_M = multi_stage(varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin == 0
    return % maybe do some default bahviour
elseif nargin == 1 % only params given 
    params = varargin{1};
elseif nargin == 2 % params and map
    params = varargin{1};
    setup     = varargin{2};
elseif nargin == 3 % params, map and where to store 
    params    = varargin{1};
    setup     = varargin{2};
    rslts     = varargin{3};
else % maybe error
    return
end

mu  = params.sim.mu;
M = 5;
N_tot = nan;
N_stg = params.sim.N_stg;
K_stg = floor(N_stg/mu);
K_tot = nan;

% Energy parameters
E_total = setup.total_energy; % [J]
E_m     = E_total;
E_min   = calc_real_energy(K_stg, params.sim.V_max*ones(1,N_stg), params);

%% Simulation Setup
% Basestation
S_b = setup.base_station_pos;
S_s = S_b;
% Communication user
S_c = setup.comm_user_pos;
% Sensing target
S_t = setup.sense_target_pos;

% estimate the target randomly 
S_target_est = setup.est_sense_target; %S_t + [100; -100];
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

% Inital trajectory
[S_init, V_init] = init_trajectory(S_s, S_mid, N_stg, params);

% calc initial values of parametes
delta_square_init = sqrt(1 + norms(V_init, 2, 1).^4/(4*params.energy.v_0^4)) - norms(V_init, 2 ,1).^2/(2*params.energy.v_0^2);
xi_init = delta_square_init;

% run the mth stage
[S_opt_m,E_m_used, V_m, xi_m, delta_m,CRB_vec_m,R_vec_m] = single_stage(E_m, N_stg, delta_square_init,K_stg, S_c, S_init,S_target_est,S_s,V_init,params);

D_meas(:,m) = sense_target(S_t, S_opt_m(:,mu:mu:end));

% store calculated trajectory 
S_opt_mat(:,1:size(S_opt_m,2), m) = S_opt_m;

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

rslts.R_opt_vecs = R_opt_vecs;
rslts.CRB_opt_vecs  = CRB_opt_vecs;
rslts.xi_m_vec     = xi_m_vec;
rslts.delta_m_vec = delta_m_vec;
rslts.E_used_vec = E_used_vec;
rslts.E_m_vec = E_m_vec;
rslts.E_min_vec = E_min_vec;
rslts.S_target_est_mat = S_target_est_mat;
rslts.S_init_mat = S_init_mat;                
rslts.V_m_mat = V_m_mat;                
rslts.K_tot = K_tot;                
rslts.N_tot = N_tot;                
rslts.M = M;                
                   
results_M.R_opt_vecs = R_opt_vecs;
results_M.CRB_opt_vecs  = CRB_opt_vecs;
results_M.xi_m_vec     = xi_m_vec;
results_M.delta_m_vec = delta_m_vec;
results_M.E_used_vec = E_used_vec;
results_M.E_m_vec = E_m_vec;
results_M.E_min_vec = E_min_vec;
results_M.S_target_est_mat = S_target_est_mat;
results_M.S_init_mat = S_init_mat;                
results_M.V_m_mat = V_m_mat;                
results_M.K_tot = K_tot;                
results_M.N_tot = N_tot;                
results_M.M = M;             
end

