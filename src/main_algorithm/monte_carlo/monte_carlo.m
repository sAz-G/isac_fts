function results_mc = monte_carlo(params, setup, number_mc_iterations)
%MONTO_CARLO Summary of this function goes here
%   Detailed explanation goes here
CRB_mc = zeros(1, number_mc_iterations);
Rate_mc = zeros(1, number_mc_iterations);

%% Do the main Monte-Carlo-Simulation
for mc_iter = 1:number_mc_iterations
    fprintf('   Monte-Carlo: n = %.f/%.f\n', mc_iter, number_mc_iterations);
    % generate random communication and sensing target
    setup.comm_user_pos = [params.sim.L_x * rand(1); params.sim.L_y * rand(1)];
    setup.sense_target_pos = [params.sim.L_x * rand(1); params.sim.L_y * rand(1)];

    % run the multi_stage algorithm for a random setup
    res = multi_stage(params, setup);
    
    mc_est_sense_target = res.S_target_est_mat(:,res.M);

    mc_traj = reshape(res.S_opt_mat(:, :, 1:res.M), 2,[]); % check when incomplete last stage

    mc_hover_traj = res.S_opt_mat(:, params.sim.mu:params.sim.mu:params.sim.N_stg, :); % check when incomplete last stage
    mc_hover_traj = reshape(mc_hover_traj(:, :, 1:res.M), 2, []); % check when incomplete last stage

    number_setpoints = length(mc_traj(1, :)); % check when incomplete last stage

    CRB_mc(mc_iter) = crb(mc_hover_traj, setup.sense_target_pos, params);
    Rate_mc(mc_iter) = avg_data_rate(mc_traj, setup.comm_user_pos, params, number_setpoints);
end

results_mc = struct('CRB_avg', {mean(CRB_mc)}, ...
                    'Rate_avg', {mean(Rate_mc)});
end