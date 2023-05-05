function [S_init, V_init] = init_trajectory(s_start, s_end,N, params)
%INIT_TRAJECTORY computes the initial trajectory for the iterative process
% arguments: 
%   s_start - starting point, the point at which the quad starts.
%   N       - the amount of trajectory points at the mth stage.
%   s_mid   - the middle point between the communication user and the
%       estimated sensing target (the estimated position obtained from the
%       radar).
%   V_str   - the speed of the quad. It is constant throughout the
%   initial trajectory.
%   T_f     - flight time between points.
% output:
%   S_init  - the initial trajectory.
%   V_out   - the velocity throughout the trajectory.

V_str = params.sim.V_str;
T_f   = params.sim.T_f;

S_init = zeros(2,N);

% inital velocity
V_init  = ones(2,N).*V_str.*(s_end - s_start)./norm(s_end - s_start);
% initial trajectory of the drone
S_init(:,1:N) = s_start + V_str .*(1:N).* (s_end - s_start) ./ norm(s_end - s_start)*T_f;
end

