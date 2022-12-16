function R = dt_rate(n,P,sigma_0,B)
%RN_DASH Summary of this function goes here
%   Detailed explanation goes here
R = B*log2(1 + (P*comm_dist(n))./(sigma_0.^2));
end

