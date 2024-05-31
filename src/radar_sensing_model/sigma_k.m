%------------------------------------------------------------------------
% FUNCTION NAME: sigma_k
% AUTHOR: Sharif Azem     (sAz-G on github), Markus Krantzik (mardank on github)
%
% DESCRIPTION: Calculate the std of the distance measurement
%
% INPUTS:
%   d_s     - distance to the sensing target
%   params  - sensing parameters
%
% OUTPUTS:
%   sig_k - standard deviation of the measurement 
%
% USAGE:  sig_k = sigma_k(d_s,params)
%
%------------------------------------------------------------------------

function sig_k = sigma_k(d_s,params)
    
    % sensing parameters 
    P              = params.sim.P;
    G_p            = params.sim.G_p;
    beta_0         = params.sim.beta_0;
    a              = params.sim.a;
    sigma_0        = params.sim.sigma_0;
    
    % standard deviation 
    sig_k = sqrt( (a.*sigma_0.^2)./( P.*G_p.*g_k(d_s, beta_0) ) );
end