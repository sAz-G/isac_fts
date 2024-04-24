%------------------------------------------------------------------------
% FUNCTION NAME: g_k
% AUTHORS: Sharif Azem (sAz-G on GitHub), Markus Krantzik (mardank on GitHub)
%
% DESCRIPTION: 
%   The two-way channel power gain.
%
% INPUTS:
%   d_s - Distance to sensing target.
%   beta_0 - Channel power at the reference distance d_s = 1m.
%
% OUTPUTS:
%   y - Channel power gain.
%
% USAGE:   
%   y = g_k(d_s, beta_0)
%
%------------------------------------------------------------------------

function y = g_k(d_s, beta_0)
    y = beta_0 ./ (d_s .^ 4);
end
