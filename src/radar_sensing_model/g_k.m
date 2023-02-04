function y = g_k(d_s,beta_0)
%G_K Summary of this function goes here
%   Detailed explanation goes here
% beta_0   = db2pow(-47); % channel power gain in dB at the reference distanc d_s(k) = 1m [W]

y = beta_0./(d_s.^4);
end

