%------------------------------------------------------------------------
% SCRIPT NAME: parameters
% AUTHORS: Sharif Azem (sAz-G on GitHub), Markus Krantzik (mardank on GitHub)
%
% DESCRIPTION: 
%   a script to set predefined parameters
%
%------------------------------------------------------------------------


% predefined parameters for the algorithm
clear 
clc
rng default;
params = struct('sim', ...
                   {struct( ... % parameters from table 1 in the paper and other few params
                           'alpha_0'   ,{db2pow(-50)},                  ...
                           'N_0'       ,{db2pow(-170)*1e-3},            ...
                           'P'         ,{db2pow(20)*1e-3},              ...
                           'B'         ,{1e06},                         ...
                           'H'         ,{200},                          ...
                           'V_str'     ,{25},                         ...
                           'L_x'       ,{1500},                         ...
                           'L_y'       ,{1500},                         ...
                           'T_f'       ,{1.5},                          ...
                           'T_h'       ,{1},                            ...
                           'beta_0'    ,{db2pow(-47)},                  ...
                           'V_max'     ,{30},                           ...
                           'mu'        ,{5},                            ...
                           'eta'       ,{0.5},                            ...
                           'iter'      ,{20},                           ...
                           'a'         ,{10},                           ...
                           'w_star'    ,{.8},                            ... %.5
                           'N_stg'     ,{25}                            ...
                          )}, ... % end of simulation parameters (sim)
                'energy', ...
                   {struct(...
                           'P_0'       , {80},      ...
                           'U_tip'     , {120},     ...
                           'D_0'       , {0.6},     ...
                           'rho'       , {1.225},   ...
                           'P_I'       , {88.6},    ...
                           'v_0'       , {4.03},    ...
                           's'         , {0.05},    ...
                           'A'         , {0.503}    ...
                           )}, ...
                 'opt_settings', ...
                            {struct(...
                            'solver'                , {'sdpt3'},          ...  % diefferent solvers - differnet reults. overall similar
                            'precision'             , {'medium'},           ...  % higher not necassarily better
                            'threshold'             , {1e-20},            ...  % lower might be worse. try out to find best.
                            'screen_out'            , {0}                ...  % screen output     
                            )}...
               ); % end params 

% additional parameters that depend on other parameters                          
params.sim.K_stg    = floor(params.sim.N_stg./params.sim.mu);
params.sim.G_p      = 0.1.*params.sim.B;
params.sim.sigma_0  = sqrt(params.sim.B.*params.sim.N_0);

setup = struct('base_station_pos', {[100; 100]},            ...
               'comm_user_pos',    {[1300; 1200]},          ...
               'est_sense_target', {[1200; 700]},           ...
               'sense_target_pos', {[200; 1300]},           ...
               'total_energy',     {35e3}                   ...
               );

