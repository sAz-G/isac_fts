function [derivatey_theta_a_out, derivatey_theta_b_out, derivatey_theta_c_out] = compute_gradient_theta_y(S_hover_in, S_target_est_in, d_s_km_in)
%COMPUTE_GRADIENT_THETA_Y Summary of this function goes here
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

    % Calculate the gradient of theta_a with respect to the y vector
    derivatey_theta_a_out = factor_CRB * x_target_est_diff.^2 .* (- 6 * y_target_est_diff) ./ (d_s_km_in.^8); % vectorized
    derivatey_theta_a_out = derivatey_theta_a_out + 8 * x_target_est_diff.^2 .* (- 4 * y_target_est_diff) ./ (d_s_km_in.^6); % vectorized

    % Calculate the gradient of theta_b with respect to the y vector
    derivatey_theta_b_part1 = (2* y_target_est_diff .* d_s_km_in.^2 - 6  * y_target_est_diff.^3) ./ (d_s_km_in.^8); % vectorized
    derivatey_theta_b_part2 = (2* y_target_est_diff .* d_s_km_in.^2 - 4  * y_target_est_diff.^3) ./ (d_s_km_in.^6); % vectorized
    derivatey_theta_b_out = factor_CRB * derivatey_theta_b_part1 + 8 * derivatey_theta_b_part2;

    % Calculate the gradient of theta_c with respect to the y vector
    derivatey_theta_c_out = factor_CRB * x_target_est_diff .* (d_s_km_in.^2 - 6 * y_target_est_diff.^2)./(d_s_km_in.^8); % vectorized
    derivatey_theta_c_out = derivatey_theta_c_out + 8 * x_target_est_diff .* (d_s_km_in.^2 - 4 * y_target_est_diff.^2)./(d_s_km_in.^6); % vectorized
end