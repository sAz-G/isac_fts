%%
clear;
clc;
close all;

if ispc
    addpath(genpath("..\..\..\src"));
elseif isunix
    addpath(genpath("../../../src"));
end 

% set rng to default

%% Simulation parameter

% load the simulation parameters 
run("call_hyperParam.m")

number_mc_iterations = 16;
L_x = params.sim.L_x;
L_y = params.sim.L_y;

energy_vec = 10e3:5e3:35e3;
CRB_over_energy = zeros(1,length(energy_vec));
Rate_over_energy = zeros(1,length(energy_vec));

counter = 1;

%% Variation of the total energy
for cur_total_energy = energy_vec
    fprintf('Variation : n = %.f/%.f, Energy: %.f\n', counter, length(energy_vec), cur_total_energy);
    % set the total available energy
    setup.total_energy = cur_total_energy;

    % call the monte_carlo function to to the monte-carlo-simulation of the current setup
    res_mc = monte_carlo(params, setup, number_mc_iterations);
    
    % save the avg. CRB
    CRB_over_energy(counter) = res_mc.CRB_avg;

    % save the avg. Rate
    Rate_over_energy(counter) = res_mc.Rate_avg;
    counter = counter + 1;
end

%% Evaluate and plot the parameters
out_path = create_output_path(fullfile('monte_carlo_variations','var_energy'));

figure
yyaxis left
semilogy(energy_vec, CRB_over_energy)
grid on

xlabel("E_{total} ")
title("MC-Simulation over E_{total} (iter=" + params.sim.iter + " eta=" + params.sim.eta + ')')
ylabel("avg. CRB")

yyaxis right
semilogy(energy_vec, Rate_over_energy)
ylabel("avg. Rate")

% save the figure
saveas(gcf, fullfile(out_path, 'MC_CRB_Rate_over_energy'), 'epsc');
saveas(gcf, fullfile(out_path, 'MC_CRB_Rate_over_energy'), 'png');

% Safe the workspace
save(fullfile(out_path, 'mc_var_energy.mat'));