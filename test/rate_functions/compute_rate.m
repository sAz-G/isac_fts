function gradient_rate_x_out = compute_rate(S_traj_in, S_c, H_in, N_in)
%COMPUTE_GRADIENT Summary of this function goes here
%   Detailed explanation goes here

    P = db2mag(20e-3); % transmit power [dB]
    N_0 = db2pow(-170e-3); % [dB/Hz]
    B = 1e6; % bandwidth [Hz]
    alpha_0 = db2pow(-50); % channel power at reference distance d_c(n) = 1m [dB]
    sigma_0 = sqrt(N_0 * B); % noise power [dB]

    d_c = compute_dc(S_traj_in, S_c, H_in, N_in);
    
    h_n = alpha_0 ./ (d_c.^2);
    gradient_rate_x_out = 1/N_in *sum(B * log2(1 + (P * h_n)/(sigma_0^2)));
end

