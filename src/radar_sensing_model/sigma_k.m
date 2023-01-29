function sig_k = sigma_k(d_s)
    N_0     = db2pow(-170) * 1e-3;  % [W/Hz]
    B       = 1e06;                 % channel bandwidth [Hz]
    sigma_0 = sqrt(B * N_0);        % unitless. Noise power at the receiver 
    P       = db2pow(20) * 1e-3;    % transmit power [dB]
    G_p      = 0.1 * B;     % signal processing gain [Hz]
    a = 10; % pre-determined constant related to the system setting

    sig_k = (a*sigma_0^2)./(P*G_p*g_k(d_s));
end
