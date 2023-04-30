%% multi_stage_script is an etnry point for the multi stage algorithm.
% it loads the parameters that are needed, which include the simulation
% setup and parameters that are given in the paper. The initial values that
% are needed for each stage are calculated and passed to the function
% optimize_m* which optimizes the trajectory at the mth stage. 
% All variables that are considered in the optimization process in each
% stage are stored in matrices. These variables include the optimization
% functions and functions and variables that are considered in the
% constraints. Also the estimated sensing target positions obtained at each
% stage is stored in a matrix. 
%
% the function multi_stage is based on this script.
%%
clear;
clc;
close all;

addpath(genpath("..\..\..\src")); % add all project folders to path  


%% Simulation parameter
run("call_hyperParam.m")  % load all predefined constans. script is in helper_functions

% call the multistage function
res = multi_stage(params, setup);

S_opt_mat = res.S_opt_mat;
s_b = setup.base_station_pos; % base station position 
s_c = setup.comm_user_pos; % position of the communication user 
s_t = setup.sense_target_pos;        % real position of the sensing target
S_target_est_mat = [setup.est_sense_target, res.S_target_est_mat];

% last stage here
[~,~,~,~,~,~,f1] = plot_map(S_opt_mat, s_b, s_t, S_target_est_mat, s_c,params);  % plot map 
f2 = figure(2);
%plot(res.J')
plot(res.J(~isnan(res.J(:,1)),~isnan(res.J(1,:)))');
ylabel('$\eta\widetilde{CRB}_\mathrm{Taylor}^m - (1 - \eta)\overline{R}^m_{\mathrm{Taylor}}$', 'Interpreter','latex');
xlabel('Iterations')
title('Oscillation at different stages')
grid on
legend("m=1","m=2","m=3","m=4","m=5","m=6")
legend('Location','southeast')
sv= 0;

if sv
    workspace_name = "ener" + setup.total_energy + "_iter" + params.sim.iter + "_omega" + params.sim.w_star + "_eta" + params.sim.eta; 
    workspace_name = strrep(workspace_name, '.','');
    out_path = create_output_path(fullfile('oscillation',workspace_name));
    saveas(f1, fullfile(out_path, 'map'), 'png');
    saveas(f1, fullfile(out_path, 'map'), 'epsc');
    saveas(f2, fullfile(out_path, 'osc'), 'png');
    saveas(f2, fullfile(out_path, 'osc'), 'epsc');
    save(fullfile(out_path, workspace_name + ".mat"));
end