clear;
close all;

%% Simulation parameter

% import the constant parameters
run("call_hyperParam.m")
eta = 1; 
mu  = 5; %

%% Simulation Setup
% Basestation
S_s = [100; 100];
% Communication user
S_c = [1.3e3; 1.2e3];
% Sensing target
S_t = [200; 1.3e3];

N_tot = 70;

%% Multi-stage approach for UWV trajectory design - First stage

% Number of setpoints for the drone
K_tot = floor(N_tot/mu); % Number of hover points

% Energy parameters
E_total = 40e3; % [J]
E_m = E_total;

% target estimation via grid search
S_target_est = S_t + L_x/20*randn(2,1);

% Middle point between communication user und target user
S_mid = (S_c + S_target_est)/2;

% Inital trajectory
[S_traj_init, V_init] = init_trajectory(S_s, N_tot, S_mid, V_str, T_f);
plot_map(S_traj_init, S_s, S_t, S_c);

%d_c_last = sqrt(H^2 + norms((S_traj_init - S_c), 2, 1))
d_c_last = norms([(S_traj_init - S_c); H*ones(1, N_tot)], 2, 1);
omega_c_last = (P * alpha_0) ./ d_c_last;
delta_square_last = sqrt(1 + norms(V_init, 2, 1).^4/(4*v_0^4)) - norms(V_init, 2 ,1).^2/(2 * v_0^2);

%% Optimization
cvx_begin
    variable S(2,N_tot)
    variable V(2,N_tot)
    variable delta(1,N_tot)
    variable xi(1,N_tot)
  
    %% Calculation of the CRB
    S_hover = S_traj_init(:, mu:mu:N_tot);

    CRB_const = compute_crb(S_hover, S_target_est, H, K_tot);

    S_hover_x = S_traj_init(1,mu:mu:N_tot);
    S_hover_y = S_traj_init(2,mu:mu:N_tot);
    
    %% Derivative with respect to x
    derivatex_CRB = compute_gradient_x(S_hover, S_target_est, H, K_tot);

    x_opt = S(1,mu:mu:N_tot); % optimization variable
    CRB_derivate_x_km = sum(derivatex_CRB .* (x_opt - S_hover_x)); % vectorized
    
    %% Derivative with respect to y
    derivatey_CRB = compute_gradient_y(S_hover, S_target_est, H, K_tot);

    y_opt = S(2,mu:mu:N_tot); % optimization variable
    CRB_derivate_y_km = sum(derivatey_CRB .* (y_opt - S_hover_y)); % vectorized
    
    % sum_up for total CRB
    CRB_affine = CRB_const + CRB_derivate_x_km + CRB_derivate_y_km;

    %% Average Communication Rate
    R_const = compute_rate(S_traj_init, S_c, H, N_tot);
    
    %% Derivative with respect to x
    derivatex_rate = compute_gradient_rate_x(S_traj_init, S_target_est, H, N_tot);

    x_traj_opt = S(1,:); % optimization variable
    R_derivate_x = sum(derivatex_rate .* (x_traj_opt - S_traj_init(1,:))); % vectorized
    
    %% Derivative with respect to y
    derivatey_rate = compute_gradient_rate_y(S_traj_init, S_target_est, H, N_tot);

    y_traj_opt = S(2,:); % optimization variable
    R_derivate_y = sum(derivatey_rate .* (y_traj_opt -  S_traj_init(1,:))); % vectorized

    %R_affine = sum(B/N_tot * 1/log(2) * sum(1./(1 + omega_c_last)) .* (omega_c - omega_c_last));
    R_affine = R_const + R_derivate_x + R_derivate_y;

    %% Objective function
    % minimization
    minimize(eta * CRB_affine - (1 - eta) * R_affine);

    %% Conditions
    Em_sum1 = sum(P_0 * (1 + 3 * pow_pos(norms(V,2,1), 2)/(U_tip.^2)) + 0.5 * D_0 * rho * s * A * pow_pos(norms(V, 2, 1), 3));
    Em_sum2 = sum(P_I*delta);
    Em_sum3 = K_tot* (P_0 + P_I);

    subject to
        E_total >= T_f * Em_sum1 + T_f * Em_sum2 + T_h * Em_sum3;
            
        S(:, 1) == S_s;

        V_max >= norms(V, 2, 1);
        delta >= 0;
        S >= 0;
        S >= 0;
        L_x >= S;
        L_y >= S;

        xi >= 0;
        % omega_c >= 0;

        % omega_c >= P * alpha_0/sigma_0^2 * pow_pos(inv_pos(d_c), 2);

        % H^2 + norms((S_traj_init - S_c), 2, 1) + diag((S_traj_init - S_c).' * (S - S_traj_init)).' >= pow_pos(d_c, 2);
        
        for i = 1:N_tot
            
            norm(V_init(:, i))^2 / v_0^2 + 2/v_0^2 * V_init(:, i).' * (V(:,i) - V_init(:,i)) >= 1/delta_square_last(i) - xi(i);

            delta_square_last(i) + 2 * sqrt(delta_square_last(i)) * (delta(i) - sqrt(delta_square_last(i))) >= xi(i);
        end
        
cvx_end

plot_map(S, S_s, S_t, S_c);
