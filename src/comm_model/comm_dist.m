function dist = comm_dist(s,x,h)
%COMM_DIST Computes the distance to the communication user
%   h == height of the drone
%   s == position of the communicaiton user
%   x == Position of the drone

dist = sqrt(h.^2+ norm(x-s).^2); 
end

