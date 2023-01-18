%%
clear;
clc;
close all;

%% Simulation parameter
run("call_hyperParam.m")
eta = 1; 
mu  = 5; 
%% Simulation Setup
% Basestation
S_s = [100; 100];
% Communication user
S_c = [1.3e3; 1.2e3];
% Sensing target
S_t = [200; 1.3e3];

N_tot = 25;

%% Multi-stage approach for UWV trajectory design - First stage

% Number of setpoints for the drone
K_tot = floor(N_tot/mu); % Number of hover points

% Energy parameters
E_total = 40e3; % [J]
E_m = E_total;

% target estimation via grid search
S_target_est = S_t ;

% Middle point between communication user und target user
S_mid = (S_c + S_target_est)/2;

% Inital trajectory
[S_init, V_init] = init_trajectory(S_s, N_tot, S_mid, V_str, T_f);
plot_map(S_init, S_s, S_t, S_c);

%d_c_last = sqrt(H^2 + norms((S_traj_init - S_c), 2, 1))
d_c_last = norms([(S_init - S_c); H*ones(1, N_tot)], 2, 1);
omega_c_last = (P * alpha_0) ./ d_c_last;
delta_square_last = sqrt(1 + norms(V_init, 2, 1).^4/(4*v_0^4)) - norms(V_init, 2 ,1).^2/(2 * v_0^2);


%% Optimization
w_star = 0.2;
iter = 5;

for u = 1:iter
% Calculation of the CRB
S_hover = S_init(:, mu:mu:N_tot);

CRB_const = compute_crb(S_hover, S_target_est, H, K_tot);

S_hover_x = S_init(1,mu:mu:N_tot);
S_hover_y = S_init(2,mu:mu:N_tot);

% Derivative with respect to x
derivatex_CRB = compute_gradient_x(S_hover, S_target_est, H, K_tot);

% Derivative with respect to y
derivatey_CRB = compute_gradient_y(S_hover, S_target_est, H, K_tot);

%y_opt = S(2,mu:mu:N_tot); % optimization variable

cvx_begin
    variable S(2,N_tot)
    variable V(2,N_tot)
    variable delta(1,N_tot)
    variable xi(1,N_tot)
    variable d_c(1,N_tot)
    variable omega_c(1,N_tot)

 
    %x_opt = S(1,mu:mu:N_tot); % optimization variable
    CRB_derivate_x_km = sum(derivatex_CRB .* (S(1,mu:mu:N_tot) - S_hover_x)); % vectorized
    
 
    CRB_derivate_y_km = sum(derivatey_CRB .* (S(2,mu:mu:N_tot) - S_hover_y)); % vectorized
    
    % sum_up for total CRB
    CRB_affine = CRB_const + CRB_derivate_x_km + CRB_derivate_y_km;

    %% Average Communication Rate

    R_affine = sum(B/N_tot * 1/log(2) * sum(1./(1 + omega_c_last)) .* (omega_c - omega_c_last));
    
    % sum up over omega_km

    %% Objective function
    % minimization
    minimize(eta * CRB_affine - (1 - eta) * R_affine);

    %% Conditions
    Em_sum1 = sum(P_0 * (1 + 3 * pow_pos(norms(V,2,1), 2)/(U_tip.^2)) + 0.5 * D_0 * rho * s * A * pow_pos(norms(V, 2, 1), 3));
    Em_sum2 = sum(P_I*delta);
    Em_sum3 = K_tot* (P_0 + P_I);

    subject to
        E_total >= T_f * Em_sum1 + T_f * Em_sum2 + T_h * Em_sum3

        V_max >= norms(V, 2, 1);
        delta >= 0;
        S(:,1) == S_s;
        S(1,:) >= 0;
        S(2,:) >= 0;
        L_x >= S(1,:);
        L_y >= S(2,:);

        xi >= 0;
        omega_c >= 0;
        
        % pow_pos(d_c, 2) >= P * alpha_0 * inv_pos(omega_c);
        omega_c >= P * alpha_0/sigma_0^2 * pow_pos(inv_pos(d_c), 2);

        % norms([(S - S_c .* ones(2,N_tot)); H*ones(1, N_tot)], 2, 1) >= d_c;
        H^2 + norms((S_init - S_c), 2, 1) + diag((S_init - S_c).' * (S - S_init)).' >= pow_pos(d_c, 2);
        
        for i = 1:(N_tot-1)
            norm(S(:,i+1) - S(:,i)) <= V_max*T_f;
            norm(V_init(:, i))^2 / v_0^2 + 2/v_0^2 * V_init(:, i).' * (V(:,i) - V_init(:,i)) >= 1/delta_square_last(i) - xi(i);

            delta_square_last(i) + 2 * sqrt(delta_square_last(i)) * (delta(i) - sqrt(delta_square_last(i))) >= xi(i);
        end
        
cvx_end

S_init = S_init + w_star.*(S-S_init);
end
plot_map(S, S_s, S_t, S_c);
