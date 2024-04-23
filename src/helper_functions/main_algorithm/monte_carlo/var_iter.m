%------------------------------------------------------------------------
% SCRIPT: var_iter
% AUTHORS: Sharif Azem (sAz-G on GitHub) and Markus Krantzik (mardank on GitHub)
%
% DESCRIPTION: Performs Monte Carlo simulations on iteration values.
%
% INPUTS:
%   Calls parameters script.
%
% OUTPUTS:
%   Saves workspace and plots.
%
% USAGE: Run the script.
%
%------------------------------------------------------------------------
%%
clear;
clc;
close all;

% Add path depending on the operating system
if ispc
    addpath(genpath("..\..\..\src"));
elseif isunix
    addpath(genpath("../../../src"));
end 

%% Simulation parameter

% Load the simulation parameters 
run("call_hyperParam.m")

number_mc_iterations = 35;
L_x = params.sim.L_x;
L_y = params.sim.L_y;

iter_vec        = 0:5:30;
iter_vec(1)     = 1;
CRB_over_iter   = zeros(1,length(iter_vec));
Rate_over_iter  = zeros(1,length(iter_vec));
MSE             = zeros(1,length(iter_vec));

counter = 1;

% Generate positions
start_bound = 0;
setup.comm_user_pos    = [start_bound + (params.sim.L_x - start_bound) * rand(1, number_mc_iterations); 
                           start_bound + (params.sim.L_y - start_bound) * rand(1, number_mc_iterations)];
setup.sense_target_pos = [start_bound + (params.sim.L_x - start_bound) * rand(1, number_mc_iterations); 
                           start_bound + (params.sim.L_y - start_bound) * rand(1, number_mc_iterations)];
setup.est_sense_target = [(params.sim.L_x) * rand(1, number_mc_iterations); 
                           (params.sim.L_y) * rand(1, number_mc_iterations)];

%% Variation of the iterations
for cur_iter = iter_vec
    fprintf('Variation : n = %d/%d, Iterations: %.f\n', counter, length(iter_vec), cur_iter);
    % Set the number of iterations
    params.sim.iter = cur_iter;

    % Call the monte_carlo function to do the Monte Carlo simulation of the current setup
    res_mc = monte_carlo(params, setup, number_mc_iterations);
    
    % Save the avg. CRB
    CRB_over_iter(counter) = mean(res_mc.CRB_avg);
    MSE(counter) = mean(res_mc.MSE_avg);

    % Save the avg. Rate
    Rate_over_iter(counter) = mean(res_mc.Rate_avg);
    counter = counter + 1;
end

%% Evaluate and plot the parameters
out_path = create_output_path(fullfile('monte_carlo_variations', 'var_iter'));

figure
title("MC-Simulation over iter")
plot(iter_vec, (CRB_over_iter))
grid on
xlabel("iter")
ylabel("CRB")

saveas(gcf, fullfile(out_path, 'MC_CRB_MSE_over_iter'), 'epsc');
saveas(gcf, fullfile(out_path, 'MC_CRB_MSE_over_iter'), 'png');

figure 
plot(iter_vec, Rate_over_iter)
xlabel("iter")
ylabel("avg. Rate")
grid on

% Save the figure
saveas(gcf, fullfile(out_path, 'MC_Rate_over_iter'), 'epsc');
saveas(gcf, fullfile(out_path, 'MC_Rate_over_iter'), 'png');

% Save the workspace
save(fullfile(out_path, 'mc_var_iter.mat'));
