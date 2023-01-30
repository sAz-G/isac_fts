clear;
close all;

%% Simulation parameter
eta = .5;

% import the constant parameters
run("call_hyperParam.m")

mu = 5; 

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
S_target_est = S_t + L_x/20*randn(2,1);

% Middle point between communication user und target user
S_mid = (S_c + S_target_est)/2;

% Inital trajectory
[S_traj_init, V_init] = init_trajectory(S_s, N_tot, S_mid, V_str, T_f);

S_last_iter = S_traj_init;
V_last_iter = V_init;
% d_c_last = norms([(S_last_iter - S_c); H*ones(1, N_tot)], 2, 1);
% omega_c_last = (P * alpha_0) ./ d_c_last;
delta_square_last = sqrt(1 + norms(V_last_iter, 2, 1).^4/(4*v_0^4)) - norms(V_last_iter, 2 ,1).^2/(2 * v_0^2);

counter_iterations = 0;

while counter_iterations < 10
    %% Optimization
    cvx_begin
        variable S(2,N_tot)
        variable V(2,N_tot-1)
        variable delta(1,N_tot)
        variable xi(1,N_tot)
    
        %% Calculation of the CRB
        S_hover = S_last_iter(:, mu:mu:N_tot);
    
        CRB_const = compute_crb(S_hover, S_target_est, H, K_tot);
    
        S_hover_x = S_last_iter(1,mu:mu:N_tot);
        S_hover_y = S_last_iter(2,mu:mu:N_tot);
        
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
        V == (S(:,2:N_tot) - S(:,1:(N_tot-1)))./T_f;

        subject to
            E_total >= T_f * Em_sum1 + T_f * Em_sum2 + T_h * Em_sum3
    
            V_max >= norms(V, 2, 1);
            delta >= 0;
            S >= 0;
            S >= 0;
            S(:,1) == S_s;
            L_x >= S;
            L_y >= S;
    
            xi >= 0;
            % omega_c >= 0;
            
            % omega_c >= P * alpha_0/sigma_0^2 * pow_pos(inv_pos(d_c), 2);
    
            % H^2 + norms((S_last_iter - S_c), 2, 1) + diag((S_last_iter - S_c).' * (S - S_last_iter)).' >= pow_pos(d_c, 2);
            
            for i = 1:N_tot-1  
                norm(V_last_iter(:, i))^2 / v_0^2 + 2/v_0^2 * V_last_iter(:, i).' * (V(:,i) - V_last_iter(:,i)) >= 1/delta_square_last(i) - xi(i);
    
                delta_square_last(i) + 2 * sqrt(delta_square_last(i)) * (delta(i) - sqrt(delta_square_last(i))) >= xi(i);
                
            end
            
    cvx_end
    %w_star = .5;
    S_last_iter = S;
    %S_last_iter = S_last_iter + w_star.*(S-S_last_iter);

    V_last_iter = V;
    delta_square_last = delta.^2;   

    counter_iterations = counter_iterations + 1;
end

plot_map(S, S_s, S_t,S_target_est, S_c);

