function gradient_rate_y_out = compute_gradient_rate_y(S_target_in, S_c, H_in, N_in)
%COMPUTE_GRADIENT Summary of this function goes here
%   Detailed explanation goes here

    % import the constant parameters
    run("call_hyperParam.m")
    
    d_c = compute_dc(S_target_in, S_c, H_in, N_in);
    
    y_c_diff = (S_target_in(2,:) - S_c(2));
    
    derivatex_dc = y_c_diff ./ d_c;
    derivatex_inner_log = (P * alpha_0)/(sigma_0^2) * (-2)./(d_c .^ 3) .* derivatex_dc;
    gradient_rate_y_out = B/(N_in * log(2)) * 1 ./ (1 + (P * alpha_0)/(sigma_0^2) * 1./(d_c.^2)) .* derivatex_inner_log;
end