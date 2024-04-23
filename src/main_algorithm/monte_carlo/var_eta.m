%------------------------------------------------------------------------
% SCRIPT: var_eta
% AUTHORS: Sharif Azem (sAz-G on GitHub), Markus Krantzik (mardank on GitHub)
%
% DESCRIPTION:
%   This script conducts Monte Carlo simulations to analyze the impact of varying the parameter eta
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

% Define a vector of eta values to explore
eta_vec = 1:-.2:0;

% Preallocate arrays to store results
CRB_over_eta  = zeros(1,length(eta_vec));
Rate_over_eta = zeros(1,length(eta_vec));
MSE           = zeros(1,length(eta_vec));

counter = 1;  % Counter for iteration

% Generate random positions for communication users, sensing targets, and estimated sensing targets
start_bound = 0;
setup.comm_user_pos    = [start_bound+ (L_x-start_bound)*rand(1,number_mc_iterations) ; start_bound + (L_y-start_bound)*rand(1,number_mc_iterations) ];
setup.sense_target_pos = [start_bound+ (L_x-start_bound)*rand(1,number_mc_iterations) ; start_bound + (L_y-start_bound)*rand(1,number_mc_iterations) ];
setup.est_sense_target = [ (L_x)*rand(1,number_mc_iterations); (L_y)*rand(1,number_mc_iterations) ];

%% Variation of the parameter eta
for cur_eta = eta_vec
    fprintf('Variation : n = %.f/%.f, eta: %.g\n', counter, length(eta_vec), cur_eta);

    % Set the current eta value in the simulation parameters
    params.sim.eta = cur_eta;

    % Perform Monte Carlo simulation for the current setup
    res_mc = monte_carlo(params, setup, number_mc_iterations);
    
    % Calculate and store the average CRB, MSE, and rate
    CRB_over_eta(counter) = mean(res_mc.CRB_avg);
    MSE(counter)           = mean(res_mc.MSE_avg);
    Rate_over_eta(counter) = mean(res_mc.Rate_avg);
    
    counter = counter + 1;  % Increment counter
end

%% Evaluate and plot the parameters
out_path = create_output_path(fullfile('monte_carlo_variations','var_eta'));

% Plot CRB and MSE against eta
figure
title("MC-Simulation over eta")
semilogy(eta_vec, CRB_over_eta)
grid on
xlabel("eta")
ylabel("")
hold on 
semilogy(eta_vec, MSE)
hold off
xlabel("eta")
ylabel("Estimation Error")
legend("CRB", "MSE")

% Save the plot
saveas(gcf, fullfile(out_path, 'MC_CRB_MSE_over_eta'), 'epsc');
saveas(gcf, fullfile(out_path, 'MC_CRB_MSE_over_eta'), 'png');

% Plot average rate against eta
figure 
plot(eta_vec, Rate_over_eta)
xlabel("eta")
ylabel("avg. Rate")
grid on

% Save the plot
saveas(gcf, fullfile(out_path, 'MC_Rate_over_eta'), 'epsc');
saveas(gcf, fullfile(out_path, 'MC_Rate_over_eta'), 'png');

% Save workspace variables
save(fullfile(out_path, 'mc_var_eta.mat'));
