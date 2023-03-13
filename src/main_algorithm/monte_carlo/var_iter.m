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

iter_vec = [5, 10, 25, 50, 100];
CRB_over_iter = zeros(1,length(iter_vec));
Rate_over_iter = zeros(1,length(iter_vec));

counter = 1;

%% Variation of the iterations
for cur_iter = iter_vec
    fprintf('Variation : n = %.f/%.f, Iterations: %.f\n', counter, length(iter_vec), cur_iter);
    % set the nmber of iterations
    parmas.sim.iter = cur_iter;

    % call the monte_carlo function to to the monte-carlo-simulation of the current setup
    res_mc = monte_carlo(params, setup, number_mc_iterations);
    
    % save the avg. CRB
    CRB_over_iter(counter) = res_mc.CRB_avg;

    % save the avg. Rate
    Rate_over_iter(counter) = res_mc.Rate_avg;
    counter = counter + 1;
end

%% Evaluate and plot the parameters
out_path = create_output_path(fullfile('monte_carlo_variations','var_iter'));

figure
yyaxis left
semilogy(iter_vec, CRB_over_iter)
grid on

xlabel("number of iterations")
title("MC-Simulation over iterations: (eta=" + params.sim.eta + ")")
ylabel("avg. CRB")

yyaxis right
semilogy(iter_vec, Rate_over_iter)
ylabel("avg. Rate")

% save the figure
saveas(gcf, fullfile(out_path, 'MC_CRB_Rate_over_iter'), 'epsc');
saveas(gcf, fullfile(out_path, 'MC_CRB_Rate_over_iter'), 'png');

% Safe the workspace
save(fullfile(out_path, 'MC_CRB_Rate_over_energy.mat'));