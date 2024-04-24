%------------------------------------------------------------------------
% FUNCTION NAME: grid_vectors
% AUTHORS: Sharif Azem (sAz-G on GitHub), Markus Krantzik (mardank on GitHub)
%
% DESCRIPTION: 
%   Create vectors for the grid to use grid search on.
%
% INPUTS:
%   L_x - The limit in x direction.
%   L_y - The limit in y direction.
%   r_x - Amount of grid points in x direction.
%   r_y - Amount of grid points in y direction.
%
% OUTPUTS:
%   x_t - The output vector for the x direction.
%   y_t - The output vector for the y direction.
%
% USAGE: 
%   [x_t, y_t] = grid_vectors(Lx, Ly, rx, ry)
%
%------------------------------------------------------------------------

function [x_t, y_t] = grid_vectors(Lx, Ly, rx, ry)
    x_t = linspace(0, Lx, rx);
    y_t = linspace(0, Ly, ry);
end
