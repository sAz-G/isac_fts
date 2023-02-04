function sig_k = sigma_k(d_s,params)
    
    P              = params.sim.P;
    G_p            = params.sim.G_p;
    beta_0         = params.sim.beta_0;
    a              = params.sim.a;
    sigma_0        = params.sim.sigma_0;
    
    sig_k = sqrt( (a.*sigma_0.^2)./( P.*G_p.*g_k(d_s, beta_0) ) );
end
