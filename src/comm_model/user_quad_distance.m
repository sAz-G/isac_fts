%------------------------------------------------------------------------
% FUNCTION NAME: user_quad_distance
% AUTHORS: Sharif Azem (sAz-G on GitHub) and Markus Krantzik (mardank on GitHub)
%
% DESCRIPTION: Calculates the distance to the communication user.
%
% INPUTS:
%   S   - Trajectory.
%   s_c - Position of the communication user.
%   H   - Height.
%
% OUTPUTS:
%   relative_distance - Distance vector.
%
% USAGE: relative_distance = user_quad_distance(S, s_c, H)
%-----------------------------------------------------------------------

function relative_distance = user_quad_distance(S, s_c, H)
    relative_pos_x = (S(1, :) - s_c(1));
    relative_pos_y = (S(2, :) - s_c(2));
    N = size(S, 2);
    relative_distance = norms([relative_pos_x; relative_pos_y; H * ones(1, N)], 2, 1);
end
