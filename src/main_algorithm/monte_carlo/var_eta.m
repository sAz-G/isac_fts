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

eta_vec = 10e3:5e3:35e3;
CRB_over_eta = zeros(1,length(eta_vec));
Rate_over_eta = zeros(1,length(eta_vec));

counter = 1;

%% Variation of the parameter eta
for cur_eta = eta_vec
    fprintf('Variation : n = %.f/%.f, eta: %.f\n', counter, length(eta_vec), cur_eta);
    % set the total available energy
    params.sim.eta = cur_eta;

    % call the monte_carlo function to to the monte-carlo-simulation of the current setup
    res_mc = monte_carlo(params, setup, number_mc_iterations);
    
    % save the avg. CRB
    CRB_over_eta(counter) = res_mc.CRB_avg;

    % save the avg. Rate
    Rate_over_eta(counter) = res_mc.Rate_avg;
    counter = counter + 1;
end

%% Evaluate and plot the parameters
out_path = create_output_path(fullfile('monte_carlo_variations','var_eta'));

figure
yyaxis left
semilogy(eta_vec, CRB_over_eta)
grid on

xlabel("eta")
title("MC-Simulation over eta: (iter=" + params.sim.iter +")")
ylabel("avg. CRB")

yyaxis right
semilogy(eta_vec, Rate_over_eta)
ylabel("avg. Rate")

% save the figure
saveas(gcf, fullfile(out_path, 'MC_CRB_Rate_over_eta'), 'epsc');
saveas(gcf, fullfile(out_path, 'MC_CRB_Rate_over_eta'), 'png');

% Safe the workspace
save(fullfile(out_path, 'mc_var_eta.mat'));