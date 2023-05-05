%------------------------------------------------------------------------
% FUNCTION NAME: grid_vectors
% AUTHOR: Sharif Azem     (TU-Darmstadt department 18, sAz-G on github)
%         Markus Krantzik (TU-Darmstadt department 18, mardank on github)
%
% DESCRIPTION: % create vectors for the grid to use grid search on
%
% INPUTS:
%   L_x - the limit in x direction.
%   L_y - the limit in y direction.
%   r_x - amount of grid points.
%
% OUTPUTS:
%   x_t - the output vector for the x direction
%   y_t - the output vector for the y direction
%
% USAGE:  [x_t,y_t] = grid_vectors(Lx,Ly,rx,ry)
%
%------------------------------------------------------------------------

function [x_t,y_t] = grid_vectors(Lx,Ly,rx,ry)
x_t = linspace(0,Lx,rx);
y_t = linspace(0,Ly,ry);
end

