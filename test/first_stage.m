clear;
close all;

%% Simulation parameter
nu = 0.5;

alpha_0 = db2pow(-50); % channel power at reference distance d_c(n) = 1m [dB]
beta_0 = db2pow(-47); % channel power at reference distance d_s(k) = 1m [dB]
N_0 = db2pow(-170e-3); % [dB/Hz]
B = 1e6; % bandwidth [Hz]
sigma_0 = sqrt(N_0 * B); % noise power [dB]
P = db2mag(20e-3); % transmit power [dB]
G_p = 0.1 * B; % Signal processing gain [Hz]
V_max = 30; % maximum velocity [m/s]
H = 200; % flying height [m]
T_h = 1; % hovering duration [s]
V_str = 20; % [m/s]
mu = 5; %
L_x = 1500; % dimension of the ground in x direction [m]
L_y = 1500; % dimension of the ground in y direction [m]
T_f = 1.5; % flying duration [s]

a = 10; % pre-determined constant related to the system setting

% Energy
W = 20; % aircraft weight in Newton
delta_eng = 0.012;
rho = 1.225; % air density [kg/m^3]
A = 0.503; % rotor_disc area in [m^2]
Sigma = 300; % blade angular velocity in [radians/second]
s = 0.05; % rotor solidity [m^3]
R = 0.4; % rotor radius in [m]
k = 0.1; % incremental correction factor to include power;

% P_0 = delta_eng/8 * rho * s * A * Sigma^3 * R^3;
P_0 = 80; % blade profile power [W]

% P_I = (1 + k) * (W^(3/2))/sqrt(2*rho * A);
P_I = 88.6; % induced power in hovering status [W]
U_tip = 120; % tip speed of the rotor [m/s]
v_0 = 4.03; % mean rotor induced velocity in forward flying [m/s]
D_0 = 0.6; % Fuselage drag ratio


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
plot_trajectory(S_traj_init, S_s, S_t, S_c)

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
    variable d_c(1,N_tot)
    variable omega_c(1,N_tot)

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

    R_affine = sum(B/N_tot * 1/log(2) * sum(1./(1 + omega_c_last)) .* (omega_c - omega_c_last));
    
    % sum up over omega_km

    %% Objective function
    % minimization
    minimize(nu * CRB_affine - (1 - nu) * R_affine);

    %% Conditions
    Em_sum1 = sum(P_0 * (1 + 3 * pow_pos(norms(V,2,1), 2)/(U_tip.^2)) + 0.5 * D_0 * rho * s * A * pow_pos(norms(V, 2, 1), 3));
    Em_sum2 = sum(P_I*delta);
    Em_sum3 = K_tot* (P_0 + P_I);

    subject to
        E_total >= T_f * Em_sum1 + T_f * Em_sum2 + T_h * Em_sum3

        V_max >= norms(V, 2, 1);
        delta >= 0;
        S >= 0;
        S >= 0;
        L_x >= S;
        L_y >= S;

        xi >= 0;
        omega_c >= 0;
        
        % pow_pos(d_c, 2) >= P * alpha_0 * inv_pos(omega_c);
        omega_c >= P * alpha_0/sigma_0^2 * pow_pos(inv_pos(d_c), 2);

        % norms([(S - S_c .* ones(2,N_tot)); H*ones(1, N_tot)], 2, 1) >= d_c;
        H^2 + norms((S_traj_init - S_c), 2, 1) + diag((S_traj_init - S_c).' * (S - S_traj_init)).' >= pow_pos(d_c, 2);
        
        for i = 1:N_tot
            
            norm(V_init(:, i))^2 / v_0^2 + 2/v_0^2 * V_init(:, i).' * (V(:,i) - V_init(:,i)) >= 1/delta_square_last(i) - xi(i);

            delta_square_last(i) + 2 * sqrt(delta_square_last(i)) * (delta(i) - sqrt(delta_square_last(i))) >= xi(i);
        end
        
cvx_end
