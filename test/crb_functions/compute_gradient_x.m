function [gradient_x_out] = compute_gradient_x(S_hover_in, S_target_est_in, H_in, K_tot_in)
%COMPUTE_GRADIENT Summary of this function goes here
%   Detailed explanation goes here

    P = db2mag(20e-3); % transmit power [dB]
    N_0 = db2pow(-170e-3); % [dB/Hz]
    B = 1e6; % bandwidth [Hz]
    G_p = 0.1 * B; % Signal processing gain [Hz]
    beta_0 = db2pow(-47); % channel power at reference distance d_s(k) = 1m [dB]
    a = 10; % pre-determined constant related to the system setting
    sigma_0 = sqrt(N_0 * B); % noise power [dB]

    [theta_a, theta_b, theta_c] = compute_thetas(S_hover_in, S_target_est_in, H_in, K_tot_in);
    CRB_denominator = theta_a * theta_b - theta_c^2;
    
    factor_CRB = (P*G_p*beta_0)/(a*sigma_0^2);
    
    x_target_est_diff = (S_hover_in(1,:) - S_target_est_in(1));
    y_target_est_diff = (S_hover_in(2,:) - S_target_est_in(2));

    % d_s_km = sqrt(H^2 + x_target_est_diff.^2 + y_target_est_diff.^2);
    d_s_km = norms([x_target_est_diff; y_target_est_diff; H_in*ones(1, K_tot_in)], 2, 1);

    % Derivative with respect to x
    derivatex_theta_a_part1 = (2* x_target_est_diff .* d_s_km.^2 - 6  * x_target_est_diff.^3)/(d_s_km.^8); % vectorized
    derivatex_theta_a_part2 = (2* x_target_est_diff .* d_s_km.^2 - 4  * x_target_est_diff.^3)/(d_s_km.^6); % vectorized

    derivatex_theta_a = factor_CRB * derivatex_theta_a_part1 + 8 * derivatex_theta_a_part2;

    derivatex_theta_b = factor_CRB * y_target_est_diff.^2 .* (- 6 * x_target_est_diff)./(d_s_km.^8); % vectorized
    derivatex_theta_b = derivatex_theta_b + 8 * y_target_est_diff.^2 .* (- 4 * x_target_est_diff)./(d_s_km.^6); % vectorized

    derivatex_theta_c = factor_CRB * y_target_est_diff .* (d_s_km.^2 - 6 * x_target_est_diff.^2)./(d_s_km.^8); % vectorized
    derivatex_theta_c = derivatex_theta_c + 8 * y_target_est_diff .* (d_s_km.^2 - 4 * x_target_est_diff.^2)./(d_s_km.^6); % vectorized
    
    derivatex_CRB_denominator = derivatex_theta_a .* theta_b + theta_a .* derivatex_theta_b - 2 * theta_c .* derivatex_theta_c;  % vectorized

    derivatex_CRB_part1 = (derivatex_theta_a .* CRB_denominator - theta_a .* derivatex_CRB_denominator)./(CRB_denominator.^2); % vectorized
    derivatex_CRB_part2 = (derivatex_theta_b .* CRB_denominator - theta_b .* derivatex_CRB_denominator)./(CRB_denominator.^2); % vectorized
    
    gradient_x_out = (derivatex_CRB_part1 + derivatex_CRB_part2); % vectorized
    
end

