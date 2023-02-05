function Drate_Dy = compute_gradient_rate_y(s_t, s_c, H, N,params)
%COMPUTE_GRADIENT Summary of this function goes here
%   Detailed explanation goes here

    % import the constant parameters
    P              = params.sim.P;
    sigma_0        = params.sim.sigma_0;
    alpha_0        = params.sim.alpha_0;
%     B              = params.sim.B;
    
    d_c = compute_dc(s_t, s_c, H, N);
    
    y_c_diff = (s_t(2,:) - s_c(2));
    
    derivatex_dc = y_c_diff ./ d_c;
    derivatex_inner_log = (P * alpha_0)/(sigma_0^2) * (-2)./(d_c .^ 3) .* derivatex_dc;
    Drate_Dy = 1/(N * log(2)) * 1 ./ (1 + (P * alpha_0)/(sigma_0^2) * 1./(d_c.^2)) .* derivatex_inner_log;
end