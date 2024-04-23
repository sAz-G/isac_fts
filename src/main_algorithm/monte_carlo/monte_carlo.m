%------------------------------------------------------------------------
% FUNCTION NAME: monte_carlo
% AUTHORS: Sharif Azem (sAz-G on GitHub) and Markus Krantzik (mardank on GitHub)
%
% DESCRIPTION: Performs Monte Carlo simulations.
%
% INPUTS:
% params               - Parameters.
% setup                - Setup parameters.
% number_mc_iterations - Number of Monte Carlo iterations.
%
% OUTPUTS:
% results_mc - Results of the Monte Carlo simulations.
%
% USAGE: results_mc = monte_carlo(params, setup, number_mc_iterations)
%------------------------------------------------------------------------

function results_mc = monte_carlo(params, setup, number_mc_iterations)

MSE_mc = zeros(1, number_mc_iterations);
CRB_mc = zeros(1, number_mc_iterations);
Rate_mc = zeros(1, number_mc_iterations);

tic
parfor mc_iter = 1:number_mc_iterations
    fprintf('   Monte-Carlo: n = %.f/%.f\n', mc_iter, number_mc_iterations);
    [mc_traj, mc_hover_traj, last_pos_est] = run_simulation(params, setup);
    
    % Estimation error
    MSE_mc(mc_iter) = compute_MSE(last_pos_est, setup.sense_target_pos(:, mc_iter));
    
    % CRB
    CRB_mc(mc_iter) = compute_CRB(mc_hover_traj, setup.sense_target_pos, params);
    
    % Data rate
    Rate_mc(mc_iter) = compute_data_rate(mc_traj, setup.comm_user_pos, params);
end
mc_runtime = toc;

results_mc = struct('CRB_avg',  {CRB_mc},  ...
                    'Rate_avg', {Rate_mc}, ...
                    'MSE_avg',  {MSE_mc},  ...
                    'Runtime_avg', {mc_runtime/number_mc_iterations});

end

function [mc_traj, mc_hover_traj, last_pos_est] = run_simulation(params, setup)
    % Generate random communication and sensing target
    setup.comm_user_pos = setup.comm_user_pos(:, mc_iter);
    setup.sense_target_pos = setup.sense_target_pos(:, mc_iter);
    setup.est_sense_target = setup.est_sense_target(:, mc_iter);

    % Run the multi-stage algorithm for a random setup
    res = multi_stage(params, setup);
    
    % Extract trajectory and hover points
    mc_traj = reshape(res.S_opt_mat(:, :, 1:res.M), 2, []);
    mc_hover_traj = res.S_opt_mat(:, params.sim.mu:params.sim.mu:params.sim.N_stg, :);
    mc_hover_traj = reshape(mc_hover_traj(:, :, 1:res.M), 2, []);
    
    % Extract last estimated position
    last_pos_est = res.S_target_est_mat(:, ~isnan(res.S_target_est_mat(1,:)));
    last_pos_est = last_pos_est(:, end);
end

function MSE = compute_MSE(last_pos_est, sense_target_pos)
    MSE = mean(norms(last_pos_est - sense_target_pos, 2, 1).^2);
end

function CRB_val = compute_CRB(mc_hover_traj, sense_target_pos, params)
    CRB_val = crb(mc_hover_traj, sense_target_pos, params);
end

function data_rate = compute_data_rate(mc_traj, comm_user_pos, params)
    data_rate = avg_data_rate(mc_traj, comm_user_pos, params, size(mc_traj, 2));
end
