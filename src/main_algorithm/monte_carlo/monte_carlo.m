function results_mc = monte_carlo(params, setup, number_mc_iterations)
%MONTO_CARLO Summary of this function goes here
%   Detailed explanation goes here
MSE_mc = zeros(1, number_mc_iterations);
CRB_mc = zeros(1, number_mc_iterations);
Rate_mc = zeros(1, number_mc_iterations);

p = gcp('nocreate');

tic
%% Do the main Monte-Carlo-Simulation
parfor mc_iter = 1:number_mc_iterations
    setup_local  = setup;
    params_local = params;
    fprintf('   Monte-Carlo: n = %.f/%.f\n', mc_iter, number_mc_iterations);

    % generate random communication and sensing target
    setup_local.comm_user_pos = setup.comm_user_pos(:,mc_iter);
    setup_local.sense_target_pos = setup.sense_target_pos(:,mc_iter);
    setup_local.est_sense_target = setup.est_sense_target(:,mc_iter);

    % run the multi_stage algorithm for a random setup
    
    res = multi_stage(params_local, setup_local);
    
    mc_est_sense_target = res.S_target_est_mat(:,res.M);

    mc_traj = reshape(res.S_opt_mat(:, :, 1:res.M), 2,[]); % check when incomplete last stage

    mc_hover_traj = res.S_opt_mat(:, params_local.sim.mu:params_local.sim.mu:params_local.sim.N_stg, :); % check when incomplete last stage
    mc_hover_traj = reshape(mc_hover_traj(:, :, 1:res.M), 2, []); % check when incomplete last stage

    number_setpoints = length(mc_traj(1, :)); % check when incomplete last stage
    
    MSE_mc(mc_iter)  = mean(norms(res.S_target_est_mat(:, ~isnan(res.S_target_est_mat(1,:)))-setup.sense_target_pos(:,number_mc_iterations),2,1).^2);
    CRB_mc(mc_iter)  = crb(mc_hover_traj, setup_local.sense_target_pos, params);
    Rate_mc(mc_iter) = avg_data_rate(mc_traj, setup_local.comm_user_pos, params, size(mc_traj,2));
end
mc_runtime = toc;

results_mc = struct('CRB_avg',  {mean(CRB_mc)}, ...
                    'Rate_avg', {mean(Rate_mc)}, ...
                    'MSE_avg', {mean(MSE_mc)}, ...
                    'Runtime_avg', {mc_runtime/number_mc_iterations});
end