%------------------------------------------------------------------------
% FUNCTION NAME: multi_stage
% AUTHORS: Sharif Azem (sAz-G on GitHub) and Markus Krantzik (mardank on GitHub)
%
% DESCRIPTION:  
%   This function performs the multistage algorithm for trajectory optimization
%   based on the parameters provided. It calculates initial values for each
%   stage and optimizes the trajectory at each stage using the function
%   optimize_m. The optimization results, including trajectory, energy usage,
%   velocity, and other variables, are stored in matrices.
%
% INPUTS:
%   params - Predefined parameters.
%   setup  - Additional setup parameters.
%
% OUTPUTS:
%   results_M - Results of the multistage algorithm stored in a structure.
%
% USAGE: 
%   results_M = multi_stage(params, setup)
%
%-----------------------------------------------------------------------

function results_M = multi_stage(params, setup)

% Simulation parameters
mu      = params.sim.mu;      % Measurement step
M       = 7;                  % Number of stages
N_stg   = params.sim.N_stg;   % Number of flight points in each stage
K_stg   = params.sim.K_stg;   % Number of hovering points in each stage

% Simulation setup
s_b = setup.base_station_pos;      % Base station position 
s_c = setup.comm_user_pos;          % Communication user position 
s_t = setup.sense_target_pos;       % Real position of the sensing target
s_target_est = setup.est_sense_target; % Estimated target position

E_total = setup.total_energy;   % Total energy
E_m = E_total;                  % Energy at the mth stage
E_min = 7e03;                   % Minimum energy

% Variables to save
m = 1;                          % Initial stage
D_meas = nan(K_stg, M);         % Measurement matrix
S_opt_mat = nan(2, N_stg, M);   % Matrix of optimal trajectory for each stage
S_init_mat = nan(2, N_stg, M);  % Matrix of initial trajectories
S_target_est_mat = nan(2, M);   % Matrix of target estimations at each stage
E_min_vec = nan(M, 1);          % Minimum energy needed at each stage
E_m_vec = nan(M, 1);            % Available energy at each stage
E_used_vec = nan(M, 1);         % Energy used at each stage
V_m_mat = nan(2, N_stg, M);     % Optimal velocity at each stage
xi_m_vec = nan(N_stg, M);       % Optimal xi at each stage
delta_m_vec = nan(N_stg, M);    % Optimal delta at each stage
R_opt_vec = nan(M, 1);          % Optimal data rate solutions
CRB_opt_vec = nan(M, 1);        % Optimal CRB solutions
R_iter_m = nan(M, params.sim.iter); % Data rate iterations
CRB_iter_m = nan(M, params.sim.iter); % CRB iterations
J = nan(M, params.sim.iter);    % Objective function values

% Optimization
epsilon = params.sim.eta;       % Weight for calculating the initial trajectory

S_total_m = [];                 % Total trajectory

while E_min < E_m
    % Initial trajectory
    s_end = s_target_est * epsilon + s_c * (1 - epsilon);
    S_init = init_trajectory(s_b, s_end, N_stg, params);
    
    S_opt_mat(:, :, m) = S_init;
    S_total_m = reshape(S_opt_mat(:, :, 1:m), 2, N_stg * m);
    
    hover_idxs = get_S_hover(params, 1, m, S_total_m);
    
    [S_opt_m, V_m, xi_m, delta_m, CRB_vec_m, R_vec_m, CRB_m, R_m, J_m] = ...
        optimize_m(E_m, s_c, S_total_m(:, hover_idxs), S_total_m, s_target_est, s_b, params);
    
    S_opt_mat(:, :, m) = S_opt_m(:, :, end);
    D_meas(:, m) = sense_target(s_t, S_opt_mat(:, mu:mu:end, m), params);
    
    s_target_est = estimate_target(S_opt_mat(:, mu:mu:end, 1:m), D_meas(1:K_stg * m), params, 'random_gridsearch');
    
    E_m_used = calc_real_energy(S_opt_mat(:, :, m), s_b, params);
    E_m = E_m - E_m_used;
    
    E_used_vec(m) = E_m_used;
    E_m_vec(m) = E_m;
    E_min_vec(m) = E_min;
    S_init_mat(:, 1:N_stg, m) = S_init;
    S_target_est_mat(:, m) = s_target_est;
    V_m_mat(:, :, m) = V_m(:, :, end);
    delta_m_vec(:, m) = delta_m(:, end);
    xi_m_vec(:, m) = xi_m(:, end);
    CRB_opt_vec(m) = CRB_vec_m;
    R_opt_vec(m) = R_vec_m;
    CRB_iter_m(m, :) = CRB_m';
    R_iter_m(m, :) = R_m';
    J(m, :) = J_m';
    
    m = m + 1;
end

M = m - 1;
N_tot = size(S_opt_mat, 2) * size(S_opt_mat, 3); % Total number of points 
K_tot = floor(N_tot / mu);                       % Total number of hover points 

% Results structure
results_M = struct('E_used_vec', E_used_vec, 'E_m_vec', E_m_vec, 'E_min_vec', E_min_vec, ...
    'S_opt_mat', S_opt_mat, 'S_init_mat', S_init_mat, 'S_target
