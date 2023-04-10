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
plot_map(S_opt_mat, s_b, s_t, S_target_est_mat, s_c,params);  % plot map 

