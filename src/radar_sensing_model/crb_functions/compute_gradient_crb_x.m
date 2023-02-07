function gradient_x_out = compute_gradient_crb_x(S_hov, S_total_hover, s_target_est, H, K_stg, params)
%COMPUTE_GRADIENT Summary of this function goes here
%   Detailed explanation goes here
    
    % Calculate the theta_a, theta_b and theta_c
    num_total_hover_p = length(S_total_hover(1, :));
    [theta_a, theta_b, theta_c] = compute_thetas(S_total_hover, s_target_est, H, num_total_hover_p, params);
    CRB_denominator = theta_a * theta_b - theta_c^2;
    
    d_s_km = compute_ds(S_hov, s_target_est, H, K_stg);

    %% Derivative with respect to x
    % function call to calculate the gradients of thetas with respect to the x vector
    [derivatex_theta_a, derivatex_theta_b, derivatex_theta_c] = compute_gradient_theta_x(S_hov, s_target_est, d_s_km, params);

    derivatex_CRB_denominator = derivatex_theta_a .* theta_b + theta_a .* derivatex_theta_b - 2 * theta_c .* derivatex_theta_c;  % vectorized

    derivatex_CRB_part1 = (derivatex_theta_a .* CRB_denominator - theta_a .* derivatex_CRB_denominator)./(CRB_denominator.^2); % vectorized
    derivatex_CRB_part2 = (derivatex_theta_b .* CRB_denominator - theta_b .* derivatex_CRB_denominator)./(CRB_denominator.^2); % vectorized
    
    gradient_x_out = (derivatex_CRB_part1 + derivatex_CRB_part2); % vectorized
end

