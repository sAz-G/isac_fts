%------------------------------------------------------------------------
% FUNCTION NAME: calc_velocity
% AUTHOR: Sharif Azem     (TU-Darmstadt department 18, sAz-G on github)
%         Markus Krantzik (TU-Darmstadt department 18, mardank on github)
%
% DESCRIPTION:  calculates the velocity matrix
%
% INPUTS:
%  S - trajectory 
%  s_s - last point of the previous stage
%  params - predefined parameters
%
% OUTPUTS:
%       V - velocity matrix 
%
% USAGE: V = calc_velocity(S, s_s, params)
%-----------------------------------------------------------------------

function V = calc_velocity(S, s_s, params)
V = diff([s_s, S], 1, 2)./params.sim.T_f;
end

