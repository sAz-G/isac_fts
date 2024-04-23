%------------------------------------------------------------------------
% FUNCTION NAME: monte_carlo
% AUTHORS: Sharif Azem (sAz-G on GitHub) and Markus Krantzik (mardank on GitHub)
%
% DESCRIPTION: Performs Monte Carlo simulations.
%
% INPUTS:
%   params               - Parameters.
%   setup                - Setup parameters.
%   number_mc_iterations - Number of Monte Carlo iterations.
%
% OUTPUTS:
%   results_mc - Results of the Monte Carlo simulations.
%
% USAGE: results_mc = monte_carlo(params, setup, number_mc_iterations)
%
%------------------------------------------------------------------------

function results_mc = monte_carlo(params, setup, number_mc_iterations)

MSE_mc = zeros(1, number_mc_iterations);
CRB_mc = zeros(1, number_mc_iterations);
Rate_mc = zeros(1, number_mc_iterations);

p = gcp('nocreate');

tic
%% Do the main Monte-Carlo-Simulation
parfor mc_iter = 1:number_mc_iterations % Run on multiple cores 
    setup_local  = setup;
    params_local = params;
    fprintf('   Monte-Carlo: n = %d/%d\n', mc_iter, number_mc_iterations);

    % Generate random communication and sensing target
    setup_local.comm_user_pos = setup.comm_user_pos(:, mc_iter);
    setup_local.sense_target_pos = setup.sense_target_pos(:, mc_iter);
    setup_local.est_sense_target = setup.est_sense_target(:, mc_iter);

    % Run the multi_stage algorithm for a random setup
    res = multi_stage(params_local, setup_local);
    
    % Trajectory 
    mc_traj = reshape(res.S_opt_mat(:, :, 1:res.M), 2, []); % Check when incomplete last stage
    
    % Hover points
    mc_hover_traj = res.S_opt_mat(:, params_local.sim.mu:params_local.sim.mu:params_local.sim.N_stg, :); % Check when incomplete last stage
    mc_hover_traj = reshape(mc_hover_traj(:, :, 1:res.M), 2, []); % Check when incomplete last stage
    
    % Last estimated position
    last_pos_est = res.S_target_est_mat(:, ~isnan(res.S_target_est_mat(1, :)));
    last_pos_est = last_pos_est(:, end);
    
    % Estimation error
    MSE_mc(mc_iter)  = mean(norms(last_pos_est - setup.sense_target_pos(:, mc_iter), 2, 1).^2);
    CRB_mc(mc_iter)  = crb(mc_hover_traj, setup_local.sense_target_pos, params);
    
    % Data rate
    Rate_mc(mc_iter) = avg_data_rate(mc_traj, setup_local.comm_user_pos, params, size(mc_traj, 2));
end
mc_runtime = toc;

results_mc = struct('CRB_avg',  {CRB_mc},  ...
                    'Rate_avg', {Rate_mc}, ...
                    'MSE_avg',  {MSE_mc},  ...
                    'Runtime_avg', {mc_runtime / number_mc_iterations});
end
