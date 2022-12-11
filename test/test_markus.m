clear;
close all;

%% Simulation parameter
nu = 0.5;

alpha_0 = -50; % [dB]
beta_0 = -47; % [dB]
N_0 = -170e-3; % [dB/Hz]
B = 1e6; % [Hz]
sigma_0 = N_0*beta_0; % [dB]
P = 20e-3; % [dB]
G_p = 0.1 * B; % [Hz]
V_max = 30; % [m/s]
H = 200; % [m]
T_h = 1; % [s]
V_str = 20; % [m/s]
mu = 5; %
L_x = 1500; % [m]
L_y = 1500; % [m]
T_f = 1.5; % [s]



%% Simulation Setup
% Basestation
X_s = [100; 100];
% Communication user
X_c = [1.3e3; 1.2e3];
% Sensing target
X_t = [200; 1.3e3];
% Middle point between communication user und target user
X_mid = (X_c + X_t)/2;

N_m = 10;

% Inital trajectory
S_init = init_trajectory(X_s, N_m, X_mid, V_str);
plot_trajectory(S_init, X_s, X_t, X_c)

% Variables for evaluation
S_out = S_init;
xt_out = null(1);
yt_out = null(1);

%% Multi-stage approach for UWV trajectory design
m = 1; % Iteration index

% Number of setpoints for the drone
N_stg = N_m; % Number of setpoints in one stage
K_stg = floor(N_stg/mu); % Number of hover points

% Energy parameters
E_total = 40e3; % [J]
E_m = E_total;
 
while (E_m > Emin)
    dm_s_est = zeros(K_m, 1);

    % Obtain distance estimate
    for j = 1:K_m
        dm_s_est(j) = ...
    end
        
    % target estimation via grid search
    [xt_m_est, yt_m_est] = ...;
    xt_out = [xt_out; xt_m_est];
    yt_out = [yt_out; yt_m_est];

    % Iteration step
    m = m + 1;
    
    % Calculate Em
    used_energy = ...
    E_m = E_m - used_energy;

    % Inital trajectory
    S_init = init_trajectory(X_s, N_m, X_mid, V_str);
    plot_trajectory(S_init, X_s, X_t, X_c)
    
    % Optain UAV trajectory
    
    % Initialization
    ...

    % Optimization
    % cvx_begin
    %     variable S(n)
    %     variable x(n)
    %     variable y(n)
    %     variable V(n)
    %     variable delta(n)
    %     variable xi(n)
    %     
    %     subject to
    %         Em >= T_f * ...
    %         Xi >= 0;
    %         for i = 1:N_m
    %             abs(V_last(i))^2 / v_0^2 + 2/v_0 * V_last.'* (V(i) - V_last(i))  >= 1/delta_last(i)^2 - Xi(i);
    %             delta_last(i)^2 + 2 * delta_last(i) * (delta(i) - delta_last) >= Xi(i)
    %         end
    %         
    % cvx_end

    % Save the opitmized variables
    S_out = [S_out; S];
end

%%% Last stage

% Parameter
N_M = N - N_m * m;
K_M = floor(N_M/mu);

% Optimization
% cvx_begin
%     variable S(n)
%     variable x(n)
%     variable y(n)
%     variable V(n)
%     variable delta(n)
%     variable xi(n)
%     
%     subject to
%         Em >= T_f * ...
%         Xi >= 0;
%         for i = 1:N_m
%             abs(V_last(i))^2 / v_0^2 + 2/v_0 * V_last.'* (V(i) - V_last(i))  >= 1/delta_last(i)^2 - Xi(i);
%             delta_last(i)^2 + 2 * delta_last(i) * (delta(i) - delta_last) >= Xi(i)
%         end
%         
% cvx_end

% Save the opitmized variables
S_out = [S_out; S];

% Obtain distance estimate
for j = 1:K_m
    dm_s_est(j) = ...
end

% target estimation via grid search
[xt_m_est, yt_m_est] = ...

%% Plot UAV trajectory


