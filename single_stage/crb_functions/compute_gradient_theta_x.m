function [derivatex_theta_a_out, derivatex_theta_b_out, derivatex_theta_c_out] = compute_gradient_theta_x(S_hover_in, S_target_est_in, d_s_km_in)
%COMPUTE_GRADIENT_THETA_X Summary of this function goes here
%   Detailed explanation goes here
    
    % import the constant parameters
    run("call_hyperParam.m")
    
    factor_CRB = (P*G_p*beta_0)/(a*sigma_0^2);
    
    x_target_est_diff = (S_hover_in(1,:) - S_target_est_in(1));
    y_target_est_diff = (S_hover_in(2,:) - S_target_est_in(2));

    % Calculate the gradient of theta_a with respect to the x vector
    derivatex_theta_a_part1 = (2* x_target_est_diff .* d_s_km_in.^2 - 6  * x_target_est_diff.^3) ./ (d_s_km_in.^8); % vectorized
    derivatex_theta_a_part2 = (2* x_target_est_diff .* d_s_km_in.^2 - 4  * x_target_est_diff.^3) ./ (d_s_km_in.^6); % vectorized
    derivatex_theta_a_out = factor_CRB * derivatex_theta_a_part1 + 8 * derivatex_theta_a_part2;

    % Calculate the gradient of theta_b with respect to the x vector
    derivatex_theta_b_out = factor_CRB * y_target_est_diff.^2 .* (- 6 * x_target_est_diff) ./ (d_s_km_in.^8); % vectorized
    derivatex_theta_b_out = derivatex_theta_b_out + 8 * y_target_est_diff.^2 .* (- 4 * x_target_est_diff) ./ (d_s_km_in.^6); % vectorized

    % Calculate the gradient of theta_c with respect to the x vector
    derivatex_theta_c_out = factor_CRB * y_target_est_diff .* (d_s_km_in.^2 - 6 * x_target_est_diff.^2) ./ (d_s_km_in.^8); % vectorized
    derivatex_theta_c_out = derivatex_theta_c_out + 8 * y_target_est_diff .* (d_s_km_in.^2 - 4 * x_target_est_diff.^2) ./ (d_s_km_in.^6); % vectorized
    
end

