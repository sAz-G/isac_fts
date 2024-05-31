%------------------------------------------------------------------------
% FUNCTION NAME: get_min
% AUTHOR: Sharif Azem (sAz-G on github), Markus Krantzik (mardank on github)
%
% DESCRIPTION: get_min runs the function to find the minimum based on the
% grid search approach
%
% INPUTS:
%   D_est   - Distance vector
%   x_t     - x positions
%   y_t     - y positions
%   x_jt    - x hovering positions
%   y_jt    - y hovering positions
%   params  - predefined parameters 
%
% OUTPUTS:
%   idx - matrix id of the x position 
%   idy - matrix id of the y position
%   min_val - minimum value 
%
% USAGE:  [idx, idy,min_val] = get_min(D_est,x_t,y_t,x_jt,y_jt,params)
%
%------------------------------------------------------------------------
function [idx, idy,min_val] = get_min(D_est,x_t,y_t,x_jt,y_jt,params)
    [idx, idy,min_val] = get_min_gridsearch(D_est,x_t,y_t,x_jt,y_jt,params);
end

