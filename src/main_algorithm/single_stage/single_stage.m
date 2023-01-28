function [S_m, E_used_m, V_m, xi_m, delta_m] = single_stage(H,eta,E_m,N_stg,delta_square_last,w_star,K_stg,iter,mu,S_c,S_init, S_t_est,S_s, V_init)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% parameters 
L_x   = 1500;
L_y   = 1500;
V_max = 30;

T_f = 1.5;
T_h = 1;

s = 0.05; 
A = 0.503;
rho = 1.225;
D_0 = 0.6;
U_tip = 120;
P_0  = 80;
P_I  = 88.6;

v_0 = 4.03;

B = 1e06;
R_max = 1;% 1e06*log(1 + (db2pow(20)*1e-3*db2pow(-50))./(sqrt(1e06*db2pow(-170) * 1e-3)*H.^2));

for u = 1:iter
% Calculation of the CRB
S_hover = S_init(:, mu:mu:N_stg);
S_hover_x = S_init(1,mu:mu:N_stg);
S_hover_y = S_init(2,mu:mu:N_stg);

% Derivative with respect to x
derivatex_CRB = compute_gradient_x(S_hover, S_t_est, H, K_stg);

% Derivative with respect to y
derivatey_CRB = compute_gradient_y(S_hover, S_t_est, H, K_stg);

% Derivative with respect to x
derivatex_rate = compute_gradient_rate_x(S_init, S_c, H, N_stg);

% Derivative with respect to y
derivatey_rate = compute_gradient_rate_y(S_init, S_c, H, N_stg);

%y_opt = S(2,mu:mu:N_tot); % optimization variable

cvx_begin
    cvx_solver mosek
    variable S(2,N_stg)
    variable V(2,N_stg)
    variable delta(1,N_stg)
    variable xi(1,N_stg)

 
    %x_opt = S(1,mu:mu:N_tot); % optimization variable

    CRB_derivate_x_km = sum(derivatex_CRB .* (S(1,mu:mu:N_stg) - S_hover_x)); % vectorized
    CRB_derivate_y_km = sum(derivatey_CRB .* (S(2,mu:mu:N_stg) - S_hover_y)); % vectorized
    
    % sum_up for total CRB
    CRB_affine =  CRB_derivate_x_km + CRB_derivate_y_km;

    % Average Communication Rate
    R_derivate_x = sum(derivatex_rate .* (S(1,:) - S_init(1,:))); % vectorized
    
    R_derivate_y = sum(derivatey_rate .* (S(2,:) -  S_init(2,:))); % vectorized

    R_affine =  R_derivate_x + R_derivate_y;
 

    %% Objective function
    % minimization
    minimize(eta * CRB_affine - (1 - eta) * R_affine./R_max);

    %% Conditions
    E_used = calc_energy(K_stg,P_I,P_0,V,U_tip,D_0,rho,s, A,delta,T_f,T_h);
    
    subject to
        E_m >= E_used;%T_f * Em_sum1 + T_f * Em_sum2 + T_h * Em_sum3;
        V_max >= norms(V, 2, 1);
        delta >= 0;
        S(:,1) == S_s;
        S(1,:) >= 0;
        S(2,:) >= 0;
        L_x >= S(1,:);
        L_y >= S(2,:);
        xi >= 0;
        
        for k = 1:N_stg-1
            %if k == 1
            %     V(k) == (S(:,k)-S_s)./T_f;
            %else
            %     V(:,k) == (S(:,k)-S(:,k-1))./T_f;
            %end
            norm(S(:,k+1) - S(:,k)) <= V_max*T_f;
           
        end
        
        for i = 1:N_stg
            
            norm(V_init(:, i))^2 / v_0^2 + 2/v_0^2 * V_init(:, i).' * (V(:,i) - V_init(:,i)) >= 1/delta_square_last(i) - xi(i);

            delta_square_last(i) + 2 * sqrt(delta_square_last(i)) * (delta(i) - sqrt(delta_square_last(i))) >=  xi(i);
        end
        
cvx_end

S_init = S_init + w_star.*(S-S_init);
%S_init = S;
V_init = V;
delta_square_last = delta.^2;
end

S_m = S_init;
E_used_m = E_used;
V_m = V;
xi_m = xi;
delta_m = delta;
end

