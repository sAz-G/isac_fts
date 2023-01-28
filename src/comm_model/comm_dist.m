function dist = comm_dist(s,x,h)
%COMM_DIST Summary of this function goes here
%   Detailed explanation goes here
dist = sqrt(h.^2+ norm(x-s).^2); 
end

