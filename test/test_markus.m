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

    %% Optimization
    cvx_begin
        variable S(2,n)
        variable x(n)
        variable y(n)
        variable V(n)
        variable delta(n)
        variable xi(n)
        variable d_c(n)
        variable omega_c(n)

        %% Calculation of the CRB
        factor_CRB = (P*G_p*beta_0)/(a*sigma_0^2);

        x_km = 
        x_km_last =
        x_km = 
        x_km_last =

        x_last_diff = (x_km - x_km_last);
        y_last_diff = (y_km - y_km_last);

        d_s_km = sqrt(H^2 + x_last_diff^2 + y_last_diff^2);
        
        x_t_est = 
        y_t_est = 
        
        d_s_vec = sqrt(H^2 + x_last_diff.^2 + y_last_diff.^2);
        
        x_diff_vec = S(1, :) - x_t_est;
        y_diff_vec = S(2, :) - y_t_est;

        sum_theta_a = sum_theta_a + sum(factor_CRB * x_diff_vec.^2 ./ (d_s_vec.^6) + 8 * x_diff_vec.^2 ./ (d_s_vec.^6) );
        sum_theta_b = sum_theta_b + sum(factor_CRB * y_diff_vec.^2 ./ (d_s_vec.^6) + 8 * y_diff_vec.^2 ./ (d_s_vec.^6) );
        sum_theta_c = sum_theta_c + sum(factor_CRB * x_diff_vec.*y_diff_vec ./ (d_s_vec.^6) + 8 * x_diff_vec.*y_diff_vec ./ (d_s_vec.^6));

        theta_a = sum_theta_a;
        theta_b = sum_theta_b;
        theta_c = sum_theta_c;

        CRB_denominator = theta_a * theta_b - theta_c^2;

        derivate_theta_a_part1 = (2* x_last_diff * d_s_km - 6 * x_last_diff^2 * x_last_diff^3)/(d_s_km^8);
        derivate_theta_a_part2 = (2* x_last_diff * d_s_km - 4 * x_last_diff^2 * x_last_diff^3)/(d_s_km^6);

        derivate_theta_a = factor_CRB * derivate_theta_a_part1 + 8 * derivate_theta_a_part2;

        derivate_theta_b = factor_CRB * y_last_diff^2 * (6 * x_last_diff)/(d_s_km^8);
        derivate_theta_b = derivate_theta_b + 8 * y_last_diff^2 * (4 * x_last_diff)/(d_s_km^6);

        derivate_theta_c = factor_CRB * y_last_diff * (d_s_km^2 - 6 * x_last_diff^2)/(d_s_km^8);
        derivate_theta_c = derivate_theta_c + 8 * y_last_diff * (d_s_km^2 - 4 * x_last_diff^2)/(d_s_km^6);
        
        derivate_CRB_denominator = derivate_theta_a * theta_b + theta_a * diff_theta_b - 2 * diff_theta_c;

        derivate_CRB_part1 = (derivate_theta_a * CRB_denominator - theta_a * derivate_CRB_denominator)/(CRB_denominator^2);
        derivate_CRB_part2 = (diff_theta_b * CRB_denominator - theta_b * derivate_CRB_denominator)/(CRB_denominator^2);

        CRB_derivate_x_km = derivate_CRB_part1 + derivate_CRB_part2;
        
        % same for y_km
        CRB_derivate_y_km = ...

        % sum_up for total CRB
        CRB_affine = CRB_derivate_x_km + CRB_derivate_y_km;

        %% Average Communication Rate
        sum_Nm = sum_Nm + Nm;

        R_affine = B/sum_Nm * 1/log(2) * 1./(1+omega_c);
        
        % sum up over omega_km

        %% Objective function
        % minimization
        minimize(nu * CRB_affine - (1 - nu) * R_affine);

        %% Conditions
        Em_sum1 = sum(P_0 * (1 + 3* abs(V)/(U_tip.^2)) + 0.5 * D_0*rho*s*A*abs(V).^3);
        Em_sum2 = sum(P_I*delta);
        Em_sum3 = K_m * (P_0 + P_I);

        subject to
            Em >= T_f * Em_sum1 + T_F * Em_sum2 + T_h * Em_sum3

            V_max >= abs(V);
            delta(i) >= 0;
            S >= 0;
            S >= 0;
            L_x >= S;
            L_y >= S;

            xi >= 0;
            
            d_c >= (P * alpha_0) / omega_c;

            for i = 1:N_m
                
                abs(V_last(i))^2 / v_0^2 + 2/v_0 * V_last.' * (V(i) - V_last(i)) >= 1/delta_last(i)^2 - Xi(i);

                delta_last(i)^2 + 2 * delta_last(i) * (delta(i) - delta_last) >= Xi(i);

                sqrt(H^2 + norm(S(:,i) - X_c)^2) >= d_c(i);
            end
            
    cvx_end

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


