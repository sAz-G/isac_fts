function [S_m, V_m, xi_m, delta_m, CRB_opt, R_opt] = single_stage_debug(E_m,N_stg,delta_square_last,K_stg,s_c,S_init_in, s_t_est,S_s, V_init, m, params)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% parameters 
L_x   = params.sim.L_x;
L_y   = params.sim.L_y;
V_max = params.sim.V_max;

v_0 = params.energy.v_0;

B = params.sim.B;
R_max = 1; %1e06*log(1 + (db2pow(20)*1e-3*db2pow(-50))./(sqrt(1e06*db2pow(-170) * 1e-3)*H.^2));
T_f = params.sim.T_f;

iter    = params.sim.iter;
mu      = params.sim.mu;
H       = params.sim.H;
eta     = params.sim.eta;
w_star  = params.sim.w_star;

CRB_opt = nan(1,iter);
R_opt   = nan(1,iter);

%% Debug parameter 

debug_S = zeros(2, N_stg, iter);
debug_V = zeros(2, N_stg, iter);
debug_delta = zeros(1, N_stg, iter);
debug_xi = zeros(1, N_stg, iter);

S_init = [S_s, S_init_in];

for u = 1:iter
% Calculation of the CRB
S_hover     = S_init(:, mu+1:mu:N_stg+1);
S_hover_x   = S_hover(1,:);
S_hover_y   = S_hover(2,:);

% Derivative with respect to x
derivatex_CRB = compute_gradient_crb_x(S_hover, s_t_est, H, K_stg, params);

% Derivative with respect to y
derivatey_CRB = compute_gradient_crb_y(S_hover, s_t_est, H, K_stg, params);

% Derivative with respect to x
derivatex_rate = compute_gradient_rate_x(S_init, s_c, H, N_stg+1,params);

% Derivative with respect to y
derivatey_rate = compute_gradient_rate_y(S_init, s_c, H, N_stg+1, params);

cvx_begin
    cvx_solver mosek
    cvx_precision high
    
    variable S(2,N_stg+1)
    variable V(2,N_stg)
    variable delta(1,N_stg)
    variable xi(1,N_stg)
    
%     expression V(2,N_stg)
% 
%     for k = 1:N_stg
%         if k == 1
%              V(:,k) = (S(:,k) - S_s)./T_f;
%         else
%              V(:,k) = (S(:,k) - S(:,k-1))./T_f;
%         end
%     end

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
    minimize(eta * CRB_affine - (1 - eta) * R_affine);

    %% Conditions
    E_used_constraint = calc_constraint_energy(K_stg,V,delta, params);
    
    subject to+
        E_m >= E_used_constraint; %T_f * Em_sum1 + T_f * Em_sum2 + T_h * Em_sum3;
        
        V_max >= norms(V, 2, 1);
        delta >= 0;
        
        S(:,1) == S_s;
        % V_max >= (S(:, 1) - S_s)./T_f;

        S(1,:) >= 0;
        S(2,:) >= 0;

        L_x >= S(1,:);
        L_y >= S(2,:);
        xi >= 0;
        %R_affine <= R_max;
        if u == 1
            xi == delta_square_last;
        end
        for k = 1:N_stg
            %if k == 1
            %     V(k) == (S(:,k)-S_s)./T_f;
            %else
            %     V(:,k) == (S(:,k)-S(:,k-1))./T_f;
            %end
            norm(S(:,k+1) - S(:,k)) <= V_max*params.sim.T_f;
           
        end

        for i = 1:N_stg
            
            norm(V_init(:, i))^2 / v_0^2 + 2/v_0^2 * V_init(:, i).' * (V(:,i) - V_init(:,i)) >= 1/delta_square_last(i) - xi(i);

            delta_square_last(i) + 2 * sqrt(delta_square_last(i)) * (delta(i) - sqrt(delta_square_last(i))) >=  xi(i);
        end
        
cvx_end

debug_S(:, :, u) = S(:,2:end);
debug_V(:, :, u) = V;
debug_delta(:, :, u) = delta;
debug_xi(:, :, u) = xi;

S_init = S_init + w_star.*(S-S_init);
% S_init = S;

V_init = V;
delta_square_last = delta.^2;

R_opt(u)    = R_affine;
CRB_opt(u)  = CRB_affine;
end

fig_conv = plot_convergence(debug_S, debug_V, debug_delta, debug_xi, S_s, m, params);

S_m = S_init(:, 2:end);
V_m = V;
xi_m = xi;
delta_m = delta;
end

