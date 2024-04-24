%------------------------------------------------------------------------
% FUNCTION NAME: sigma_k
% AUTHORS: Sharif Azem (sAz-G on GitHub), Markus Krantzik (mardank on GitHub)
%
% DESCRIPTION: 
%   Estimate the distance to the target based on the sensing model.
%
% INPUTS:
%   s_t - Sensing target position.
%   s_q - Hovering positions of the drone.
%   params - Constant parameters.
%
% OUTPUTS:
%   d_hat - Estimated distance to the target.
%
% USAGE:   
%   d_hat = sense_target(s_t, S_q, params);
%
%------------------------------------------------------------------------

function d_hat = sense_target(s_t, S_q, params)
    d_s = norms([s_t - S_q; ones(1, size(S_q, 2)) * params.sim.H], 2);
    d_hat = d_s + sigma_k(d_s, params) .* randn(1, length(d_s));
end
