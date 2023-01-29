function [derivatey_theta_a, derivatey_theta_b, derivatey_theta_c] = compute_gradient_theta_y(S_hover, s_target_est, ds, params)
%COMPUTE_GRADIENT_THETA_Y Summary of this function goes here
%   Detailed explanation goes here
   
    % import the constant parameters
      
    P              = params.sim.P;
    G_p            = params.sim.G_p;
    beta_0         = params.sim.beta_0;
    a              = params.sim.a;
    sigma_0        = params.sim.sigma_0;
    factor_CRB = (P*G_p*beta_0)/(a*sigma_0^2);

    x_target_est_diff = (S_hover(1,:) - s_target_est(1));
    y_target_est_diff = (S_hover(2,:) - s_target_est(2));

    % Calculate the gradient of theta_a with respect to the y vector
    derivatey_theta_a = factor_CRB * x_target_est_diff.^2 .* (- 6 * y_target_est_diff) ./ (ds.^8); % vectorized
    derivatey_theta_a = derivatey_theta_a + 8 * x_target_est_diff.^2 .* (- 4 * y_target_est_diff) ./ (ds.^6); % vectorized

    % Calculate the gradient of theta_b with respect to the y vector
    derivatey_theta_b_part1 = (2* y_target_est_diff .* ds.^2 - 6  * y_target_est_diff.^3) ./ (ds.^8); % vectorized
    derivatey_theta_b_part2 = (2* y_target_est_diff .* ds.^2 - 4  * y_target_est_diff.^3) ./ (ds.^6); % vectorized
    derivatey_theta_b = factor_CRB * derivatey_theta_b_part1 + 8 * derivatey_theta_b_part2;

    % Calculate the gradient of theta_c with respect to the y vector
    derivatey_theta_c = factor_CRB * x_target_est_diff .* (ds.^2 - 6 * y_target_est_diff.^2)./(ds.^8); % vectorized
    derivatey_theta_c = derivatey_theta_c + 8 * x_target_est_diff .* (ds.^2 - 4 * y_target_est_diff.^2)./(ds.^6); % vectorized
end