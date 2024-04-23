%------------------------------------------------------------------------
% SCRIPT: var_eta
% AUTHORS: Sharif Azem (sAz-G on GitHub) and Markus Krantzik (mardank on GitHub)
%
% DESCRIPTION: Performs Monte Carlo simulations on eta values.
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

% Set rng to default

%% Simulation parameter 

% Load the simulation parameters 
run("call_hyperParam.m")

number_mc_iterations = 35;
L_x = params.sim.L_x;
L_y = params.sim.L_y;

eta_vec = 1:-.2:0;
CRB_over_eta  = zeros(1,length(eta_vec));
Rate_over_eta = zeros(1,length(eta_vec));
MSE           = zeros(1,length(eta_vec));

counter = 1;

% Generate positions
start_bound = 0;
setup.comm_user_pos    = [start_bound + (params.sim.L_x - start_bound) * rand(1, number_mc_iterations); 
                           start_bound + (params.sim.L_y - start_bound) * rand(1, number_mc_iterations)];
setup.sense_target_pos = [start_bound + (params.sim.L_x - start_bound) * rand(1, number_mc_iterations); 
                           start_bound + (params.sim.L_y - start_bound) * rand(1, number_mc_iterations)];
setup.est_sense_target = [(params.sim.L_x) * rand(1, number_mc_iterations); 
                           (params.sim.L_y) * rand(1, number_mc_iterations)];

%% Variation of the parameter eta
for cur_eta = eta_vec
    fprintf('Variation : n = %d/%d, eta: %.g\n', counter, length(eta_vec), cur_eta);
    % Set the total available energy
    params.sim.eta = cur_eta;

    % Call the monte_carlo function to do the Monte Carlo simulation of the current setup
    res_mc = monte_carlo(params, setup, number_mc_iterations);
    
    % Save the avg. CRB
    CRB_over_eta(counter) = mean(res_mc.CRB_avg);
    MSE(counter)             = mean(res_mc.MSE_avg);
    % Save the avg. Rate
    Rate_over_eta(counter) = mean(res_mc.Rate_avg);
    counter = counter + 1;
end

%% Evaluate and plot the parameters
out_path = create_output_path(fullfile('monte_carlo_variations', 'var_eta'));

figure
title("MC-Simulation over eta")
semilogy(eta_vec, (CRB_over_eta))
grid on
xlabel("eta")
ylabel("")
hold on 
semilogy(eta_vec, (MSE))
hold off
xlabel("eta")
ylabel("Estimation Error")
legend("CRB", "MSE")

saveas(gcf, fullfile(out_path, 'MC_CRB_MSE_over_eta'), 'epsc');
saveas(gcf, fullfile(out_path, 'MC_CRB_MSE_over_eta'), 'png');

figure 
plot(eta_vec, Rate_over_eta)
xlabel("eta")
ylabel("avg. Rate")
grid on

% Save the figure
saveas(gcf, fullfile(out_path, 'MC_Rate_over_eta'), 'epsc');
saveas(gcf, fullfile(out_path, 'MC_Rate_over_eta'), 'png');

% Save the workspace
save(fullfile(out_path, 'mc_var_eta.mat'));
