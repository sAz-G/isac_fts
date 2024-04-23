%------------------------------------------------------------------------
% SCRIPT: var_energy
% AUTHORS: Sharif Azem (sAz-G on GitHub) and Markus Krantzik (mardank on GitHub)
%
% DESCRIPTION: Performs Monte Carlo simulations on energy levels.
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

% Variation of energy 
energy_vec = 10e3:5e3:35e3;
MSE = zeros(1,length(energy_vec));
CRB_over_energy = zeros(1,length(energy_vec));
Rate_over_energy = zeros(1,length(energy_vec));

counter = 1;

% Generate positions
start_bound = 0;
setup.comm_user_pos    = [start_bound + (params.sim.L_x - start_bound) * rand(1, number_mc_iterations); 
                           start_bound + (params.sim.L_y - start_bound) * rand(1, number_mc_iterations)];
setup.sense_target_pos = [start_bound + (params.sim.L_x - start_bound) * rand(1, number_mc_iterations); 
                           start_bound + (params.sim.L_y - start_bound) * rand(1, number_mc_iterations)];
setup.est_sense_target = [(params.sim.L_x) * rand(1, number_mc_iterations); 
                           (params.sim.L_y) * rand(1, number_mc_iterations)];

%% Variation of the total energy
for cur_total_energy = energy_vec % Iterate over energy
    fprintf('Variation : n = %d/%d, Energy: %.f\n', counter, length(energy_vec), cur_total_energy);
    % Set the total available energy
    setup.total_energy = cur_total_energy;
    % Call the monte_carlo function to do the Monte Carlo simulation of the current setup
    res_mc = monte_carlo(params, setup, number_mc_iterations);
    
    % Save the avg. CRB
    CRB_over_energy(counter) = mean(res_mc.CRB_avg);
    MSE(counter) = mean(res_mc.MSE_avg);
    % Save the avg. Rate
    Rate_over_energy(counter) = mean(res_mc.Rate_avg);
    counter = counter + 1;
end

%% Evaluate and plot the parameters
out_path = create_output_path(fullfile('monte_carlo_variations', 'var_ener'));

% Plot and save results
figure
title("Energy Variation")
semilogy(energy_vec, (CRB_over_energy))
hold on
semilogy(energy_vec, (MSE))
hold off
grid on
xlabel("Energy")
ylabel("Estimation Error")
legend("CRB", "MSE")

saveas(gcf, fullfile(out_path, 'MC_CRB_MSE_over_ener'), 'epsc');
saveas(gcf, fullfile(out_path, 'MC_CRB_MSE_over_ener'), 'png');

figure 
plot(energy_vec, Rate_over_energy)
xlabel("Energy")
ylabel("Avg. Rate")
grid on

% Save the figure
saveas(gcf, fullfile(out_path, 'MC_Rate_over_ener'), 'epsc');
saveas(gcf, fullfile(out_path, 'MC_Rate_over_ener'), 'png');

% Save the workspace
save(fullfile(out_path, 'mc_var_ener.mat'));
