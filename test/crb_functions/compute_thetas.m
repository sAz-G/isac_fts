function [theta_a_out, theta_b_out, theta_c_out] = compute_thetas(S_hover_in, S_target_est_in, H_in, K_tot_in)
%COMPUTE_THETAS Summary of this function goes here
%   Detailed explanation goes here

    P = db2mag(20e-3); % transmit power [dB]
    N_0 = db2pow(-170e-3); % [dB/Hz]
    B = 1e6; % bandwidth [Hz]
    G_p = 0.1 * B; % Signal processing gain [Hz]
    beta_0 = db2pow(-47); % channel power at reference distance d_s(k) = 1m [dB]
    a = 10; % pre-determined constant related to the system setting
    sigma_0 = sqrt(N_0 * B); % noise power [dB]

    factor_CRB = (P*G_p*beta_0)/(a*sigma_0^2);
    
    x_target_est_diff = (S_hover_in(1,:) - S_target_est_in(1));
    y_target_est_diff = (S_hover_in(2,:) - S_target_est_in(2));

    % d_s_km = sqrt(H^2 + x_target_est_diff.^2 + y_target_est_diff.^2);
    d_s_km = norms([x_target_est_diff; y_target_est_diff; H_in*ones(1, K_tot_in)], 2, 1);

    theta_a_out = sum(factor_CRB * x_target_est_diff.^2 ./ (d_s_km.^6) + 8 * x_target_est_diff.^2 ./ (d_s_km.^6) );
    theta_b_out = sum(factor_CRB * x_target_est_diff.^2 ./ (d_s_km.^6) + 8 * y_target_est_diff.^2 ./ (d_s_km.^6) );
    theta_c_out = sum(factor_CRB * x_target_est_diff.*y_target_est_diff ./ (d_s_km.^6) + 8 * x_target_est_diff.*y_target_est_diff ./ (d_s_km.^6));
end

