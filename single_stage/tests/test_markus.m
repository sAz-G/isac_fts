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
S_s = [100; 100];
% Communication user
S_c = [1.3e3; 1.2e3];
% Sensing target
S_t = [200; 1.3e3];
% Middle point between communication user und target user
S_mid = (S_c + S_t)/2;

N_m = 20;

% Inital trajectory
[S_init, V_init] = init_trajectory(S_s, N_m, S_mid, V_str);
plot_trajectory(S_init, S_s, S_t, S_c)

% Variables for evaluation
S_out = S_init;

%% Multi-stage approach for UWV trajectory design - First stage
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
        dm_s_est(j) = sqrt(H^2 + x_t_est_diff_old.^2 + y_t_est_diff_old.^2) + 10 * randn(1); % Replace with correct function
    end
        
    % target estimation via grid search
    pos_target_est = S_t + 10*randn(2,1);

    % Iteration step
    m = m + 1;
    
    % Calculate Em
    E_m = E_m - used_energy(V_init, K_m);

    % Inital trajectory
    S_init = init_trajectory(S_s, N_m, S_mid, V_str);
    plot_trajectory(S_init, S_s, S_t, S_c);
    
    % Optain UAV trajectory
    sumKm_hover = null(1);
    
    % Initialization
    
    %% Optimization
    cvx_begin
        variable S(2,Nm)
        variable x(n)
        variable y(n)
        variable V(n)
        variable delta(n)
        variable xi(n)
        variable d_c(n)
        variable omega_c(n)

        %% Calculation of the CRB
        factor_CRB = (P*G_p*beta_0)/(a*sigma_0^2);
        
        x_last = S_init(1,:);
        y_last = S_init(1,:);

        x_t_est = 
        y_t_est = ...
        
        sumKm_hover = ...
        pos_old_hover = ...

        x_t_est_diff_old = pos_old_hover(1,:) - ones(1, sumKm_hover) * x_t_est;
        y_t_est_diff_old = pos_old_hover(2,:) - ones(1, sumKm_hover) * y_t_est;

        d_s_vec_old = sqrt(H^2 + x_t_est_diff_old.^2 + y_t_est_diff_old.^2);
        
        sum_theta_a = sum(factor_CRB * x_t_est_diff_old.^2 ./ (d_s_vec_old.^6) + 8 * x_t_est_diff_old.^2 ./ (d_s_vec_old.^6) );
        sum_theta_b = sum(factor_CRB * y_t_est_diff_old.^2 ./ (d_s_vec_old.^6) + 8 * y_t_est_diff_old.^2 ./ (d_s_vec_old.^6) );
        sum_theta_c = sum(factor_CRB * x_t_est_diff_old.*y_t_est_diff_old ./ (d_s_vec_old.^6) + 8 * x_t_est_diff_old.*y_t_est_diff_old ./ (d_s_vec_old.^6));

        x_t_est_diff = (x_last - x_t_est);
        y_t_est_diff = (y_last - x_t_est);

        d_s_km = sqrt(H^2 + x_t_est_diff.^2 + y_t_est_diff.^2);
        
        theta_a = sum_theta_a + sum(factor_CRB * x_t_est_diff.^2 ./ (d_s_km.^6) + 8 * x_t_est_diff.^2 ./ (d_s_km.^6) );
        theta_b = sum_theta_b + sum(factor_CRB * x_t_est_diff.^2 ./ (d_s_km.^6) + 8 * y_t_est_diff.^2 ./ (d_s_km.^6) );
        theta_c = sum_theta_c + sum(factor_CRB * x_t_est_diff.*y_t_est_diff ./ (d_s_km.^6) + 8 * x_t_est_diff.*y_t_est_diff ./ (d_s_km.^6));

        CRB_denominator = theta_a .* theta_b - theta_c.^2; % vectorized

        %% Derivative with respect to x
        derivatex_theta_a_part1 = (2* x_t_est_diff .* d_s_km.^2 - 6  * x_t_est_diff.^3)/(d_s_km.^8); % vectorized
        derivatex_theta_a_part2 = (2* x_t_est_diff .* d_s_km.^2 - 4  * x_t_est_diff.^3)/(d_s_km.^6); % vectorized

        derivatex_theta_a = factor_CRB * derivatex_theta_a_part1 + 8 * derivatex_theta_a_part2;

        derivatex_theta_b = factor_CRB * y_t_est_diff.^2 .* (- 6 * x_t_est_diff)./(d_s_km.^8); % vectorized
        derivatex_theta_b = derivatex_theta_b + 8 * y_t_est_diff.^2 .* (- 4 * x_t_est_diff)./(d_s_km.^6); % vectorized

        derivatex_theta_c = factor_CRB * y_t_est_diff .* (d_s_km.^2 - 6 * x_t_est_diff.^2)./(d_s_km.^8); % vectorized
        derivatex_theta_c = derivatex_theta_c + 8 * y_t_est_diff .* (d_s_km.^2 - 4 * x_t_est_diff.^2)./(d_s_km.^6); % vectorized
        
        derivatex_CRB_denominator = derivatex_theta_a .* theta_b + theta_a .* derivatex_theta_b - 2 * derivatex_theta_c;  % vectorized

        derivatex_CRB_part1 = (derivatex_theta_a .* CRB_denominator - theta_a .* derivatex_CRB_denominator)./(CRB_denominator.^2); % vectorized
        derivatex_CRB_part2 = (derivatex_theta_b .* CRB_denominator - theta_b .* derivatex_CRB_denominator)./(CRB_denominator.^2); % vectorized
        
        x_opt = x; % optimization variable
        CRB_derivate_x_km = sum((derivatex_CRB_part1 + derivatex_CRB_part2) .* (x_opt - x_last)); % vectorized
        
        %% Derivative with respect to x
        derivatey_theta_a = factor_CRB * x_t_est_diff.^2 .* (- 6 * y_t_est_diff)./(d_s_km.^8); % vectorized
        derivatey_theta_a = derivatey_theta_a + 8 * x_t_est_diff.^2 .* (- 4 * y_t_est_diff)./(d_s_km.^6); % vectorized

        derivatey_theta_b_part1 = (2* y_t_est_diff .* d_s_km.^2 - 6  * y_t_est_diff.^3)/(d_s_km.^8); % vectorized
        derivatey_theta_b_part2 = (2* y_t_est_diff .* d_s_km.^2 - 4  * y_t_est_diff.^3)/(d_s_km.^6); % vectorized

        derivatey_theta_b = factor_CRB * derivatey_theta_b_part1 + 8 * derivatey_theta_b_part2;

        derivatey_theta_c = factor_CRB * x_t_est_diff .* (d_s_km.^2 - 6 * y_t_est_diff.^2)./(d_s_km.^8); % vectorized
        derivatey_theta_c = derivatey_theta_c + 8 * x_t_est_diff .* (d_s_km.^2 - 4 * y_t_est_diff.^2)./(d_s_km.^6); % vectorized
        
        derivatey_CRB_denominator = derivatey_theta_a .* theta_b + theta_a .* derivatey_theta_b - 2 * derivatey_theta_c;  % vectorized

        derivatey_CRB_part1 = (derivatey_theta_a .* CRB_denominator - theta_a .* derivatey_CRB_denominator)./(CRB_denominator.^2); % vectorized
        derivatey_CRB_part2 = (derivatey_theta_b .* CRB_denominator - theta_b .* derivatey_CRB_denominator)./(CRB_denominator.^2); % vectorized
        
        y_opt = y; % optimization variable
        CRB_derivate_y_km = sum((derivatey_CRB_part1 + derivatey_CRB_part2) .* (y_opt - y_last)); % vectorized
        
        % sum_up for total CRB
        CRB_affine = CRB_derivate_x_km + CRB_derivate_y_km;

        %% Average Communication Rate
        sum_Nm = sum_Nm + Nm;

        R_affine = B/sum_Nm * 1/log(2) * sum(1./(1+omega_c));
        
        % sum up over omega_km

        %% Objective function
        % minimization
        minimize(nu * CRB_affine - (1 - nu) * R_affine);

        %% Conditions
        Em_sum1 = sum(P_0 * (1 + 3* abs(V).^2/(U_tip.^2)) + 0.5 * D_0*rho*s*A*abs(V).^3);
        Em_sum2 = sum(P_I*delta);
        Em_sum3 = K_m * (P_0 + P_I);

        subject to
            Em >= T_f * Em_sum1 + T_f * Em_sum2 + T_h * Em_sum3

            V_max >= abs(V);
            delta(i) >= 0;
            S >= 0;
            S >= 0;
            L_x >= S;
            L_y >= S;

            xi >= 0;
            
            d_c >= (P * alpha_0) / omega_c;
            H^2 + norms((S_traj_init - S_c), 2, 1) + diag((S_traj_init - S_c).' * (S_traj_opt - S_traj_init)).' >= 
            for i = 1:N_m
                
                abs(V_last(i))^2 / v_0^2 + 2/v_0 * V_last.' * (V(i) - V_last(i)) >= 1/delta_last(i)^2 - Xi(i);

                delta_last(i)^2 + 2 * delta_last(i) * (delta(i) - delta_last) >= Xi(i);

                sqrt(H^2 + norm(S(:,i) - S_c)^2) >= d_c(i);
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


