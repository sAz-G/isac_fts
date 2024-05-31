%------------------------------------------------------------------------
% SCRIPT: var_iter
% AUTHOR: Sharif Azem (sAz-G on github), Markus Krantzik (mardank on github)
%
% DESCRIPTION: perform monte_carlo sims on iteration values
%
% INPUTS:
% calls parameters script
%
% OUTPUTS:
% saves workspace and plots
% USAGE: run script
%------------------------------------------------------------------------
%%
clear;
clc;
close all;

if ispc
    addpath(genpath("..\..\..\src"));
elseif isunix
    addpath(genpath("../../../src"));
end 

%% Simulation parameter

% load the simulation parameters 
run("call_parameters.m")

number_mc_iterations = 2;
L_x = params.sim.L_x;
L_y = params.sim.L_y;

iter_vec        = 0:5:10;
iter_vec(1)     = 1;
CRB_over_iter   = zeros(1,length(iter_vec));
Rate_over_iter  = zeros(1,length(iter_vec));
MSE             = zeros(1,length(iter_vec));

counter = 1;

% generate positions
start_bound = 0;
setup.comm_user_pos    = [start_bound+ (params.sim.L_x-start_bound)*rand(1,number_mc_iterations) ; start_bound + (params.sim.L_y-start_bound)*rand(1,number_mc_iterations) ];
setup.sense_target_pos = [start_bound+ (params.sim.L_x-start_bound)*rand(1,number_mc_iterations) ; start_bound + (params.sim.L_y-start_bound)*rand(1,number_mc_iterations) ];
setup.est_sense_target = [ (params.sim.L_x)*rand(1,number_mc_iterations); (params.sim.L_y)*rand(1,number_mc_iterations) ];

%% Variation of the iterations
for cur_iter = iter_vec
    fprintf('Variation : n = %.f/%.f, Iterations: %.f\n', counter, length(iter_vec), cur_iter);
    % set the nmber of iterations
    params.sim.iter = cur_iter;

    % call the monte_carlo function to to the monte-carlo-simulation of the current setup
    res_mc = monte_carlo(params, setup, number_mc_iterations);
    
    % save the avg. CRB
    CRB_over_iter(counter) =   mean(res_mc.CRB_avg);
    MSE(counter)             = mean(res_mc.MSE_avg);

    % save the avg. Rate
    Rate_over_iter(counter) = mean(res_mc.Rate_avg);
    counter = counter + 1;
end

%% Evaluate and plot the parameters
out_path = create_output_path(fullfile('monte_carlo_variations','var_iter'));

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

% save the figure
saveas(gcf, fullfile(out_path, 'MC_Rate_over_iter'), 'epsc');
saveas(gcf, fullfile(out_path, 'MC_Rate_over_iter'), 'png');

% Safe the workspace
save(fullfile(out_path, 'mc_var_iter.mat'));