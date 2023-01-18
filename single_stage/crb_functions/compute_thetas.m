function [theta_a_out, theta_b_out, theta_c_out] = compute_thetas(S_hover_in, S_target_est_in, H_in, K_tot_in)
%COMPUTE_THETAS Summary of this function goes here
%   Detailed explanation goes here

    % import the constant parameters
    run("call_hyperParam.m")
    
    factor_CRB = (P*G_p*beta_0)/(a*sigma_0^2);
    
    % d_s_km = sqrt(H^2 + x_target_est_diff.^2 + y_target_est_diff.^2);
    % d_s_km = norms([x_target_est_diff; y_target_est_diff; H_in*ones(1, K_tot_in)], 2, 1);
    d_s_km = compute_ds(S_hover_in, S_target_est_in, H_in, K_tot_in);

    x_target_est_diff = (S_hover_in(1,:) - S_target_est_in(1));
    y_target_est_diff = (S_hover_in(2,:) - S_target_est_in(2));

    theta_a_out = sum(factor_CRB * x_target_est_diff.^2 ./ (d_s_km.^6) + 8 * x_target_est_diff.^2 ./ (d_s_km.^4) );
    theta_b_out = sum(factor_CRB * y_target_est_diff.^2 ./ (d_s_km.^6) + 8 * y_target_est_diff.^2 ./ (d_s_km.^4) );
    theta_c_out = sum(factor_CRB * x_target_est_diff.*y_target_est_diff ./ (d_s_km.^6) + 8 * x_target_est_diff.*y_target_est_diff ./ (d_s_km.^4));
end

