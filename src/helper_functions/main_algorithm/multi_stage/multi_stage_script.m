%------------------------------------------------------------------------
% SCRIPT: multi_stage_script
% AUTHORS: Sharif Azem (sAz-G on GitHub) and Markus Krantzik (mardank on GitHub)
%
% DESCRIPTION: Perform the multi-stage algorithm.
%
% INPUTS: Calls the parameters script.
%
% OUTPUTS: Saves workspace and plots.
% USAGE: Run the script.
%------------------------------------------------------------------------

clear;
clc;
close all;

addpath(genpath("..\..\..\src")); % Add all project folders to the path  

%% Load parameters
run("call_parameters.m")  % Load all predefined constants.

%% Multi-stage algorithm
res = multi_stage(params, setup);

% Extract relevant variables
S_opt_mat = res.S_opt_mat;
s_b = setup.base_station_pos; % Base station position 
s_c = setup.comm_user_pos; % Position of the communication user 
s_t = setup.sense_target_pos; % Real position of the sensing target
S_target_est_mat = [setup.est_sense_target, res.S_target_est_mat];

% Plot map
[~,~,~,~,~,~,f1] = plot_map(S_opt_mat, s_b, s_t, S_target_est_mat, s_c, params);

% Plot oscillation
f2 = figure(2);
plot(res.J(~isnan(res.J(:,1)), ~isnan(res.J(1,:)))');
ylabel('$\eta\widetilde{CRB}_\mathrm{Taylor}^m - (1 - \eta)\overline{R}^m_{\mathrm{Taylor}}$', 'Interpreter','latex');
xlabel('Iterations')
title('Oscillation at different stages')
grid on
legend("m=1", "m=2", "m=3", "m=4", "m=5", "m=6", 'Location', 'southeast')

% Save workspace 
sv = 0;
if sv
    workspace_name = "ener" + setup.total_energy + "_iter" + params.sim.iter + "_omega" + params.sim.w_star + "_eta" + params.sim.eta; 
    workspace_name = strrep(workspace_name, '.', '');
    out_path = create_output_path(fullfile('oscillation',workspace_name));
    saveas(f1, fullfile(out_path, 'map'), 'png');
    saveas(f1, fullfile(out_path, 'map'), 'epsc');
    saveas(f2, fullfile(out_path, 'osc'), 'png');
    saveas(f2, fullfile(out_path, 'osc'), 'epsc');
    save(fullfile(out_path, workspace_name + ".mat"));
end
