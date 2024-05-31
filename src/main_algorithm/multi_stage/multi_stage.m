%------------------------------------------------------------------------
% FUNCTION NAME: multi_stage
% AUTHOR: Sharif Azem     (sAz-G on github), Markus Krantzik (mardank on github)
%         
%
% DESCRIPTION:  
%               passed the parameters that are needed, which include the simulation
%               setup and parameters that are given in the paper. The initial values that
%               are needed for each stage are calculated and passed to the function
%               optimize_m which optimizes the trajectory at the mth stage. 
%               All variables that are considered in the optimization process in each
%               stage are stored in a matrix. These variables include the optimization
%               functions and functions and variables that are considered in the
%               constraints. Also the estimated sensing target positions obtained at each
%               stage are stored in a matrix. 
%
% INPUTS:
%   params - predefined parameters
%   setup  - further setup parameters
%
% OUTPUTS:
%   results_M - results of the multistage algorithm
%
% USAGE: results_M = multi_stage(params, setup)
%
%-----------------------------------------------------------------------


function results_M = multi_stage(params, setup)


%% Simulation parameter
mu  = params.sim.mu;      % measurment step. Measure sensing target every mu steps
M = 7;                      
N_stg = params.sim.N_stg; % amount of flight points in every stage 
K_stg = params.sim.K_stg;  % amount of hovering points in every stage 

%% Simulation Setup
% Basestation
s_b = setup.base_station_pos;        % base station position 
s_s = s_b;                           % start position of the initial trajectory 

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
E_min = 7e03;

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
xi_m_vec         = nan(N_stg,M);                   % store the optimal xi at each stage 
delta_m_vec      = nan(N_stg,M);                   % store the optimal delta at each stage
R_opt_vec        = nan(M,1);        % store optimal data rate solutions 
CRB_opt_vec      = nan(M,1);        % store the optimalcrb soutions 
R_iter_m         = nan(M,params.sim.iter);
CRB_iter_m       = nan(M,params.sim.iter);
J                = nan(M,params.sim.iter);

%% Optimization
epsilon = params.sim.eta; % a weight for the calculation of the initial trajectory 
if params.sim.eta == 1    
    epsilon = 0.9; % do not set epsilon to 1. otherwise determinant of the crb can be 0.
end

S_total_m   = []; % trajectory points up to the mth stage included. From now called total trajectory

while E_min < E_m % break if the energy is not enough for additional N_stg points
    
    % end point between communication user und target user for the initial traj
    s_end = s_target_est*epsilon + s_c*(1-epsilon);
    
    % Inital trajectory
    S_init = init_trajectory(s_s, s_end, N_stg, params);
    
    S_opt_mat(:,:, m) = S_init;
    S_total_m         = reshape(S_opt_mat(:,:,1:m),2,N_stg.*m);
    
    % get hover points
    hover_idxs = get_S_hover(params, 1, m, S_total_m); % get hover points from the total traj 
%     plot_map(S_opt_mat, s_b, s_t, s_target_est, s_c,params);  % plot map 

    % get optimal solution
    [S_opt_m, V_m, xi_m, delta_m,CRB_vec_m,R_vec_m,CRB_m, R_m, J_m] = optimize_m(E_m, s_c, S_total_m(:,hover_idxs), S_total_m, s_target_est, s_s, params);
    
    % store calculated trajectory 
    S_opt_mat(:,:, m) = S_opt_m(:,:,end); % assign only final solution
    
    % sense target at each hover point 
    D_meas(:,m) = sense_target(s_t, S_opt_mat(:,mu:mu:end,m), params); % get new K_stg measurments
    
    % get the estimated target matrix index
    s_target_est  = estimate_target(S_opt_mat(:, mu:mu:end,1:m),D_meas(1:K_stg*m), params, 'random_gridsearch');
                                                   
    % calculate the available energy
    E_m_used = calc_real_energy(S_opt_mat(:,:, m), s_s, params);
    
    % set new current point
    s_s = S_opt_m(:,end,end);
    
    % calculate the energy 
    E_m = E_m - E_m_used; 
    
    % store variables
    E_used_vec(m)           = E_m_used;
    E_m_vec(m)              = E_m;
    E_min_vec(m)            = E_min;
    S_init_mat(:,1:N_stg, m) = S_init;
    S_target_est_mat(:,m)   = s_target_est;
    V_m_mat(:,:,m)          = V_m(:,:,end);
    delta_m_vec(:,m)        = delta_m(:,end);
    xi_m_vec(:,m)           = xi_m(:,end);
    CRB_opt_vec(m)          = CRB_vec_m;
    R_opt_vec(m)            = R_vec_m;
    CRB_iter_m(m,:)         = CRB_m';
    R_iter_m(m,:)           = R_m';
    J(m,:)                  = J_m';

    % increase the iteration variable
    m = m+1;
end

M = m-1;   

% amount of stages 
N_tot = size(S_opt_mat,2)*size(S_opt_mat,3); % total amount of points 
K_tot = floor(N_tot/mu);                     % total amount of hover points 

% struct for returning the results
results_M = struct('E_used_vec'             , {E_used_vec},         ...
                   'E_m_vec'                , {E_m_vec},            ...
                   'E_min_vec'              , {E_min_vec},          ...
                   'S_opt_mat'              , {S_opt_mat},          ...
                   'S_init_mat'             , {S_init_mat},         ...
                   'S_target_est_mat'       , {S_target_est_mat},   ...
                   'D'                      , {D_meas},             ...
                   'V_m_mat'                , {V_m_mat},            ...
                   'delta_m_vec'            , {delta_m_vec},        ...
                   'xi_m_vec'               , {xi_m_vec},           ...
                   'CRB_opt_vecs'           , {CRB_opt_vec},        ...
                   'R_opt_vecs'             , {R_opt_vec},          ...
                   'CRB_iter_m'             , {CRB_iter_m},         ...
                   'R_iter_m'               , {R_iter_m},           ...
                   'J'                      , {J},                  ...
                   'M'                      , {M},                  ...
                   'N_tot'                  , {N_tot},              ...
                   'K_tot'                  , {K_tot});