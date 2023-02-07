function [x_t,y_t] = grid_vectors(Lx,Ly,rx,ry)
% grid_vectors(Lx,Ly,rx,ry) creates two vectors in y and x direction. rx and
% ry are the resolution parameters in each dimension
% Input:
% L_x - the limit in x direction.
% L_y - the limit in y direction.
%
% Output: 
% x_t - the output vector for the x direction, which starts at x = 0 and
% ends at x = Lx. The vector has r_x elements.
% y_t - the output vector for the y direction, which starts at y = 0 and
% ends at y = Ly. The vector has r_y elements.

x_t = linspace(0,Lx,rx);
y_t = linspace(0,Ly,ry);
end

