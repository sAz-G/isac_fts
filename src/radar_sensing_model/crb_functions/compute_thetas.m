function [theta_a, theta_b, theta_c] = compute_thetas(S_hov, s_target_est, H, K_stg, params)
%COMPUTE_THETAS Summary of this function goes here
%   Detailed explanation goes here

    P              = params.sim.P;
    G_p            = params.sim.G_p;
    beta_0         = params.sim.beta_0;
    a              = params.sim.a;
    sigma_0        = params.sim.sigma_0;

    factor_CRB = (P*G_p*beta_0)/(a*sigma_0^2);
    
    d_s_km = compute_ds(S_hov, s_target_est, H, K_stg);

    x_target_est_diff = (S_hov(1,:) - s_target_est(1));
    y_target_est_diff = (S_hov(2,:) - s_target_est(2));

    theta_a = sum(factor_CRB * x_target_est_diff.^2 ./ (d_s_km.^6) + 8 * x_target_est_diff.^2 ./ (d_s_km.^4) );
    theta_b = sum(factor_CRB * y_target_est_diff.^2 ./ (d_s_km.^6) + 8 * y_target_est_diff.^2 ./ (d_s_km.^4) );
    theta_c = sum(factor_CRB * x_target_est_diff.*y_target_est_diff ./ (d_s_km.^6) + 8 * x_target_est_diff.*y_target_est_diff ./ (d_s_km.^4));
end

