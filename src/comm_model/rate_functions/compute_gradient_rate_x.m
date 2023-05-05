function gradient_rate_x_out = compute_gradient_rate_x(s_target, s_c, H, N, params)
%COMPUTE_GRADIENT_RATE_X Compute the gradient of the communication ratae with respect to x
%   s_target == position of the drone
%   s_c == position of the comm user
%   H == flying height
%   N == numer of setpoints
%   params == struct with simulation parameters
    
    % import the constant parameters
    P              = params.sim.P; % transmit power
    sigma_0        = params.sim.sigma_0; % noise power simga_0^2
    alpha_0        = params.sim.alpha_0; % channel gain at reference distance
%     B              = params.sim.B;
    
    % distance to communication user
    d_c = compute_dc(s_target, s_c, H, N); 
    
    % Compute the gradient with respect to x
    x_c_diff = (s_target(1,:) - s_c(1));
    derivatex_dc = x_c_diff ./ d_c;
    derivatex_inner_log = (P * alpha_0)/(sigma_0^2) * (-2)./(d_c .^ 3) .* derivatex_dc;
    gradient_rate_x_out = 1/(N * log(2)) * 1 ./ (1 + (P * alpha_0)/(sigma_0^2) * 1./(d_c.^2)) .* derivatex_inner_log;
end

