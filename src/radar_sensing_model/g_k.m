%------------------------------------------------------------------------
% FUNCTION NAME: g_k
% AUTHOR: Sharif Azem
%         Markus Krantzik
%
% DESCRIPTION: The two way channel power gain 
%
% INPUTS:
%   d_s     - distance to sensing target
%   beta_0  - channel power at the reference distance d_s = 1m,  
%
% OUTPUTS:
%   y - channel power gain
%
% USAGE:  y = g_k(d_s,beta_0)
%
%------------------------------------------------------------------------
function y = g_k(d_s,beta_0)
y = beta_0./(d_s.^4);
end

