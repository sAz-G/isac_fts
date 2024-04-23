%------------------------------------------------------------------------
% FUNCTION NAME: calc_velocity
% AUTHORS: Sharif Azem (sAz-G on GitHub) and Markus Krantzik (mardank on GitHub)
%
% DESCRIPTION: Calculates the velocity matrix given a trajectory.
%
% INPUTS:
%   S       - Trajectory matrix.
%   s_s     - Last point of the previous stage.
%   params  - Predefined parameters.
%
% OUTPUTS:
%   V       - Velocity matrix.
%
% USAGE: V = calc_velocity(S, s_s, params)
%-----------------------------------------------------------------------

function V = calc_velocity(S, s_s, params)
    % Calculate velocity matrix
    V = diff([s_s, S], 1, 2) ./ params.sim.T_f;
end
