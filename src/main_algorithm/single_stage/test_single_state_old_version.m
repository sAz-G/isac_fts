%%
clear;
clc;
close all;

addpath(genpath("..\..\..\src"));

%% Simulation parameter
run("call_hyperParam.m")
params.sim.eta = 0.5; 
params.sim.mu  = 5; 
eta = params.sim.eta; 
mu  = params.sim.mu; 

%% Simulation Setup
% Basestation
S_s = [100; 100];
% Communication user
S_c = [1.3e3; 1.2e3];
% Sensing target
S_t = [200; 1.3e3];

% parameters 
L_x   = params.sim.L_x;
L_y   = params.sim.L_y;
V_max = params.sim.V_max;

v_0 = params.energy.v_0;

H       = params.sim.H;

N_tot = 25;

%% Multi-stage approach for UWV trajectory design - First stage

% Number of setpoints for the drone
K_tot = floor(N_tot/mu); % Number of hover points

% Energy parameters
E_total = 35e3; % [J]
E_m = E_total;

% target estimation via grid search
s_target_est = S_t ;%+ [100; -100];

% Middle point between communication user und target user
S_mid = (S_c + s_target_est)/2;

% Inital trajectory
[S_init, V_init] = init_trajectory(S_s, S_mid, N_tot, params);
plot_map(S_init, S_s, S_t, s_target_est, S_c);

delta_square_last = sqrt(1 + norms(V_init, 2, 1).^4/(4*v_0^4)) - norms(V_init, 2 ,1).^2/(2 * v_0^2);


%% Optimization
params.sim.w_star = 0.2;
params.sim.iter = 5;
w_star  = params.sim.w_star;
iter    = params.sim.iter;

for u = 1:iter
% Calculation of the CRB
S_hover = S_init(:, mu:mu:N_tot);

CRB_const = compute_crb(S_hover, s_target_est, H, K_tot, params);
R_const = compute_rate(S_init, S_c, H, N_tot, params);

S_hover_x = S_init(1,mu:mu:N_tot);
S_hover_y = S_init(2,mu:mu:N_tot);

% Derivative with respect to x
derivatex_CRB = compute_gradient_x(S_hover, s_target_est, H, K_tot, params);

% Derivative with respect to y
derivatey_CRB = compute_gradient_y(S_hover, s_target_est, H, K_tot, params);

% Derivative with respect to x
derivatex_rate = compute_gradient_rate_x(S_init, S_c, H, N_tot, params);

% Derivative with respect to y
derivatey_rate = compute_gradient_rate_y(S_init, S_c, H, N_tot, params);

%y_opt = S(2,mu:mu:N_tot); % optimization variable

cvx_begin
    variable S(2,N_tot)
    variable V(2,N_tot)
    variable delta(1,N_tot)
    variable xi(1,N_tot)

 
    %x_opt = S(1,mu:mu:N_tot); % optimization variable
    CRB_derivate_x_km = sum(derivatex_CRB .* (S(1,mu:mu:N_tot) - S_hover_x)); % vectorized
    
    CRB_derivate_y_km = sum(derivatey_CRB .* (S(2,mu:mu:N_tot) - S_hover_y)); % vectorized
    
    % sum_up for total CRB
    CRB_affine = CRB_const + CRB_derivate_x_km + CRB_derivate_y_km;
    % CRB_affine = CRB_derivate_x_km + CRB_derivate_y_km;

    % Average Communication Rate
    %x_traj_opt = S(1,:); % optimization variable
    R_derivate_x = sum(derivatex_rate .* (S(1,:) - S_init(1,:))); % vectorized
    
    R_derivate_y = sum(derivatey_rate .* (S(2,:) -  S_init(2,:))); % vectorized

    R_affine = R_const + R_derivate_x + R_derivate_y;
 

    %% Objective function
    % minimization
    minimize(eta * CRB_affine - (1 - eta) * R_affine);
    
    P_0 = params.energy.P_0;
    P_I = params.energy.P_I;
    U_tip = params.energy.U_tip;
    D_0 = params.energy.D_0;
    rho = params.energy.rho;
    s = params.energy.s;
    A = params.energy.A;

    T_f = params.sim.T_f;
    T_h = params.sim.T_h;

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
        
    
        for u = 1:N_tot-1
        
            norm(S(:,u+1) - S(:,u)) <= V_max*T_f;
           
        end
        
        for i = 1:N_tot
            
            norm(V_init(:, i))^2 / v_0^2 + 2/v_0^2 * V_init(:, i).' * (V(:,i) - V_init(:,i)) >= 1/delta_square_last(i) - xi(i);

            delta_square_last(i) + 2 * sqrt(delta_square_last(i)) * (delta(i) - sqrt(delta_square_last(i))) >= xi(i);
        end
        
cvx_end

S_init = S_init + w_star.*(S-S_init);

% S_init = S;
V_init = V;
delta_square_last = delta.^2;
end
plot_map(S, S_s, S_t, s_target_est, S_c);
