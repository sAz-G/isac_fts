function [S_m, V_m, xi_m, delta_m, CRB_opt, R_opt] = optimize_m_debug(E_m, s_c,S_hover, S_total, delta_square_init,s_t_est,s_s, V_init, params)
%optimize_m_debug is a a function to debug the 
%   Detailed explanation goes here

% parameters 
L_x   = params.sim.L_x;
L_y   = params.sim.L_y;
V_max = params.sim.V_max;

v_0 = params.energy.v_0;

iter    = params.sim.iter;
mu      = params.sim.mu;
H       = params.sim.H;
eta     = params.sim.eta;
w_star  = params.sim.w_star;
N_stg   = params.sim.N_stg;
K_stg   = floor(N_stg/mu);

CRB_opt = nan(1,iter);
R_opt   = nan(1,iter);

%% Debug parameter 

debug_S     = zeros(2, N_stg, iter);
debug_V     = zeros(2, N_stg, iter);
debug_delta = zeros(1, N_stg, iter);
debug_xi    = zeros(1, N_stg, iter);

S_init = S_total(:,end-N_stg+1:end);

for u = 1:iter

S_hover_x                           = S_init(1, mu:mu:N_stg);
S_hover_y                           = S_init(2, mu:mu:N_stg);

S_hover_to_visit  = [S_hover_x; S_hover_y];
S_hover(:,end-K_stg+1:end)  = S_hover_to_visit; 

% Derivative with respect to x
crb_grad_x = crb_grad(S_hover, s_t_est, params, 'x'); %%%%

% Derivative with respect to y
crb_grad_y = crb_grad(S_hover, s_t_est, params, 'y'); %%%%

% Derivative with respect to x
derivatex_rate = compute_gradient_rate_x(S_init, s_c, H, N_stg,params);

% Derivative with respect to y
derivatey_rate = compute_gradient_rate_y(S_init, s_c, H, N_stg, params);

cvx_begin
    cvx_solver mosek
    cvx_precision high
    
    variable S(2,N_stg)
    variable V(2,N_stg)
    variable delta(1,N_stg)
    variable xi(1,N_stg)
    
    crb_taylor_x = sum(crb_grad_x.*(S(1,mu:mu:N_stg) - S_hover_x)); % vectorized
    crb_taylor_y = sum(crb_grad_y.*(S(2,mu:mu:N_stg) - S_hover_y)); % vectorized
    
    % sum_up for total CRB
    CRB_affine =  crb_taylor_x + crb_taylor_y;

    % Average Communication Rate
    R_derivate_x = sum(derivatex_rate .* (S(1,:) - S_init(1,:)));  % vectorized
    R_derivate_y = sum(derivatey_rate .* (S(2,:) -  S_init(2,:))); % vectorized

    R_affine =  R_derivate_x + R_derivate_y;
 

    %% Objective function
    % minimization
    minimize(eta * CRB_affine - (1 - eta) * R_affine);

    %% Conditions
    E_used_constraint = calc_constraint_energy(K_stg,V,delta, params);
    
    subject to
        E_m     >= E_used_constraint; 
        V_max   >= norms(V, 2, 1);
        delta   >= 0;
        S(1,:)  >= 0;
        S(2,:)  >= 0;
        L_x     >= S(1,:);
        L_y     >= S(2,:);
        xi      >= 0;
        R_affine <= log(1 + (db2pow(20)*1e-3*db2pow(-50))./(sqrt(1e06*db2pow(-170) * 1e-3)*H.^2));
        R_affine >= log(1 + (db2pow(20)*1e-3*db2pow(-50))./(sqrt(1e06*db2pow(-170) * 1e-3)*(H.^2+ 2*L_x.^2) ) );
 
        for k = 1:N_stg
            if k == 1
                 norm(S(:,k) - s_s)     <= V_max*params.sim.T_f;
            else
                norm(S(:,k) - S(:,k-1)) <= V_max*params.sim.T_f;
            end
        end

        for i = 1:N_stg
            norm(V_init(:, i))^2 / v_0^2 + 2/v_0^2 * V_init(:, i).' * (V(:,i) - V_init(:,i)) >= pow_pos(inv_pos(delta(i)), 2) - xi(i);
            delta_square_init(i) + 2 * sqrt(delta_square_init(i)) * (delta(i) - sqrt(delta_square_init(i))) >=  xi(i);
        end
        
cvx_end

if u < iter
    S_init = S_init + w_star.*(S-S_init);
else
    S_init = S;
end

V_init = V;
delta_square_init = delta.^2;


debug_S(:, :, u) = S_init;
debug_V(:, :, u) = V;
debug_delta(:, :, u) = delta;
debug_xi(:, :, u) = xi;

R_opt(u)    = R_affine;
CRB_opt(u)  = CRB_affine;
end

%fig_conv    = plot_convergence(debug_S, debug_V, debug_delta, debug_xi, s_s, m, params);

S_m         = S_init;
V_m         = V;
xi_m        = xi;
delta_m     = delta;
end

