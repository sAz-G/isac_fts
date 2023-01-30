function gradient_rate_x_out = compute_rate(S_traj_in, S_c, H_in, N_in, params)
%COMPUTE_GRADIENT Summary of this function goes here
%   Detailed explanation goes here

    % import the constant parameters,
    alpha_0 = params.sim.alpha_0;
    P = params.sim.P;
    sigma_0 = params.sim.sigma_0;
    B = params.sim.B;
    
    d_c = compute_dc(S_traj_in, S_c, H_in, N_in);
    
    h_n = alpha_0 ./ (d_c.^2);
    gradient_rate_x_out = 1/N_in *sum(B * log2(1 + (P * h_n)/(sigma_0^2)));
end

