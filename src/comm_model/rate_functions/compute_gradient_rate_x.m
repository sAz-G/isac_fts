function gradient_rate_x_out = compute_gradient_rate_x(s_target, s_c, H, N, params)
%COMPUTE_GRADIENT Summary of this function goes here
%   Detailed explanation goes here
    
    % import the constant parameters
  
    P              = params.sim.P;
    sigma_0        = params.sim.sigma_0;
    alpha_0        = params.sim.alpha_0;
    B              = params.sim.B;
    
    d_c = compute_dc(s_target, s_c, H, N);
    
    x_c_diff = (s_target(1,:) - s_c(1));
    
    derivatex_dc = x_c_diff ./ d_c;
    derivatex_inner_log = (P * alpha_0)/(sigma_0^2) * (-2)./(d_c .^ 3) .* derivatex_dc;
    gradient_rate_x_out = B/(N * log(2)) * 1 ./ (1 + (P * alpha_0)/(sigma_0^2) * 1./(d_c.^2)) .* derivatex_inner_log;
end

