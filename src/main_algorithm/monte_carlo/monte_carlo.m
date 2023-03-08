function results_mc = monte_carlo(params, setup, number_mc_iterations)
%MONTO_CARLO Summary of this function goes here
%   Detailed explanation goes here
CRB_mc = zeros(1, number_mc_iterations);

%% Do the main Monte-Carlo-Simulation
for mc_iter = 1:number_mc_iterations
    fprintf('   Monte-Carlo: n = %.f/%.f\n', mc_iter, number_mc_iterations);
    % generate random communication and sensing target
    setup.comm_user_pos = [params.sim.L_x * rand(1); params.sim.L_y * rand(1)];
    setup.sense_target_pos = [params.sim.L_x * rand(1); params.sim.L_y * rand(1)];
    
    % run the multi_stage algorithm for a random setup
    res = multi_stage(params, setup);
    
    CRB_mc(mc_iter) = min(res.CRB_opt_vecs(:, res.M));
end

results_mc = struct('CRB_avg', {mean(CRB_mc)});
end