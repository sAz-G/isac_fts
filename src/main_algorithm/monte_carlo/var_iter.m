%------------------------------------------------------------------------
% SCRIPT: var_iter
% AUTHORS: Sharif Azem (sAz-G on GitHub), Markus Krantzik (mardank on GitHub)
%
% DESCRIPTION:
%   This script performs Monte Carlo simulations to analyze the impact of varying the number of iterations
%   on performance metrics such as the Cramer-Rao Bound (CRB), Mean Squared Error (MSE), and average rate.
%
% INPUTS:
%   - External script: call_hyperParam.m (Loads simulation parameters)
%
% OUTPUTS:
%   - Workspace variables saved in a .mat file
%   - Plots saved as .eps and .png files
%
% USAGE:
%   Execute the script to perform Monte Carlo simulations and generate plots.
%------------------------------------------------------------------------

%%
clear;              % Clear workspace
clc;                % Clear command window
close all;          % Close all figures

% Add appropriate paths based on the operating system
if ispc
    addpath(genpath("..\..\..\src"));
elseif isunix
    addpath(genpath("../../../src"));
end 

%% Simulation Parameters 

% Load simulation parameters from an external script
run("call_hyperParam.m")

number_mc_iterations = 35;  % Number of Monte Carlo iterations
L_x = params.sim.L_x;       % Length of the simulation area along the x-axis
L_y = params.sim.L_y;       % Length of the simulation area along the y-axis

% Define a vector of iteration values to explore
iter_vec        = 0:5:30;
iter_vec(1)     = 1;

% Preallocate arrays to store results
CRB_over_iter   = zeros(1,length(iter_vec));
Rate_over_iter  = zeros(1,length(iter_vec));
MSE             = zeros(1,length(iter_vec));

counter = 1;  % Counter for iteration

% Generate random positions for communication users, sensing targets, and estimated sensing targets
start_bound = 0;
setup.comm_user_pos    = [start_bound+ (L_x-start_bound)*rand(1,number_mc_iterations) ; start_bound + (L_y-start_bound)*rand(1,number_mc_iterations) ];
setup.sense_target_pos = [start_bound+ (L_x-start_bound)*rand(1,number_mc_iterations) ; start_bound + (L_y-start_bound)*rand(1,number_mc_iterations) ];
setup.est_sense_target = [ (L_x)*rand(1,number_mc_iterations); (L_y)*rand(1,number_mc_iterations) ];

%% Variation of the Iterations
for cur_iter = iter_vec
    fprintf('Variation : n = %.f/%.f, Iterations: %.f\n', counter, length(iter_vec), cur_iter);

    % Set the number of iterations in the simulation parameters
    params.sim.iter = cur_iter;

    % Perform Monte Carlo simulation for the current setup
    res_mc = monte_carlo(params, setup, number_mc_iterations);
    
    % Calculate and store the average CRB, MSE, and rate
    CRB_over_iter(counter) = mean(res_mc.CRB_avg);
    MSE(counter)           = mean(res_mc.MSE_avg);
    Rate_over_iter(counter)= mean(res_mc.Rate_avg);
    
    counter = counter + 1;  % Increment counter
end

%% Evaluate and plot the parameters
out_path = create_output_path(fullfile('monte_carlo_variations','var_iter'));

% Plot CRB against number of iterations
figure
title("MC-Simulation over iterations")
plot(iter_vec, CRB_over_iter)
grid on
xlabel("Iterations")
ylabel("CRB")

% Save the plot
saveas(gcf, fullfile(out_path, 'MC_CRB_MSE_over_iter'), 'epsc');
saveas(gcf, fullfile(out_path, 'MC_CRB_MSE_over_iter'), 'png');

% Plot average rate against number of iterations
figure 
plot(iter_vec, Rate_over_iter)
xlabel("Iterations")
ylabel("Avg. Rate")
grid on

% Save the plot
saveas(gcf, fullfile(out_path, 'MC_Rate_over_iter'), 'epsc');
saveas(gcf, fullfile(out_path, 'MC_Rate_over_iter'), 'png');

% Save workspace variables
save(fullfile(out_path, 'mc_var_iter.mat'));
