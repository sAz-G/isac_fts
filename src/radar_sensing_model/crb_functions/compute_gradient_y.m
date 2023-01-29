function gradient_y = compute_gradient_y(S_hov, s_target_est, H, K_stg, params)
%COMPUTE_GRADIENT_Y Summary of this function goes here
%   Detailed explanation goes here

    % import the constant parameters
    run("call_hyperParam.m")
    
    [theta_a, theta_b, theta_c] = compute_thetas(S_hov, s_target_est, H, K_stg,params);
    CRB_denominator = theta_a * theta_b - theta_c^2;
    
    d_s_km = compute_ds(S_hov, s_target_est, H, K_stg);

    %% Derivative with respect to x
    use_function_compute_gradient_theta_y = true;

    if ~use_function_compute_gradient_theta_y

        P              = params.sim.P;
        G_p            = params.sim.G_p;
        beta_0         = params.sim.beta_0;
        a              = params.sim.a;
        sigma_0        = params.sim.sigma_0;
        factor_CRB = (P*G_p*beta_0)/(a*sigma_0^2);
        
        x_target_est_diff = (S_hov(1,:) - s_target_est(1));
        y_target_est_diff = (S_hov(2,:) - s_target_est(2));

        % Calculate the gradient of theta_a with respect to the y vector
        derivatey_theta_a = factor_CRB * x_target_est_diff.^2 .* (- 6 * y_target_est_diff) ./ (d_s_km.^8); % vectorized
        derivatey_theta_a = derivatey_theta_a + 8 * x_target_est_diff.^2 .* (- 4 * y_target_est_diff) ./ (d_s_km.^6); % vectorized
    
        % Calculate the gradient of theta_b with respect to the y vector
        derivatey_theta_b_part1 = (2* y_target_est_diff .* d_s_km.^2 - 6  * y_target_est_diff.^3) ./ (d_s_km.^8); % vectorized
        derivatey_theta_b_part2 = (2* y_target_est_diff .* d_s_km.^2 - 4  * y_target_est_diff.^3) ./ (d_s_km.^6); % vectorized
        derivatey_theta_b = factor_CRB * derivatey_theta_b_part1 + 8 * derivatey_theta_b_part2;
    
        % Calculate the gradient of theta_c with respect to the y vector
        derivatey_theta_c = factor_CRB * x_target_est_diff .* (d_s_km.^2 - 6 * y_target_est_diff.^2)./(d_s_km.^8); % vectorized
        derivatey_theta_c = derivatey_theta_c + 8 * x_target_est_diff .* (d_s_km.^2 - 4 * y_target_est_diff.^2)./(d_s_km.^6); % vectorized
    else
        % function call to calculate the gradients of thetas with respect to the x vector
        [derivatey_theta_a, derivatey_theta_b, derivatey_theta_c] = compute_gradient_theta_y(S_hov, s_target_est, d_s_km, params);
    end

    
    derivatey_CRB_denominator = derivatey_theta_a .* theta_b + theta_a .* derivatey_theta_b - 2 * theta_c .* derivatey_theta_c;  % vectorized

    derivatey_CRB_part1 = (derivatey_theta_a .* CRB_denominator - theta_a .* derivatey_CRB_denominator)./(CRB_denominator.^2); % vectorized
    derivatey_CRB_part2 = (derivatey_theta_b .* CRB_denominator - theta_b .* derivatey_CRB_denominator)./(CRB_denominator.^2); % vectorized
    
    gradient_y = (derivatey_CRB_part1 + derivatey_CRB_part2); % vectorized
end

