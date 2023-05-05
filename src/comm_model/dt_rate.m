function R = dt_rate(n,P,sigma_0,B)
%RT_RATE Compute the communication rate
%   n = 
R = B*log2(1 + (P*comm_dist(n))./(sigma_0.^2));
end

