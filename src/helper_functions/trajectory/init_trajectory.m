function [S_init, V_init] = init_trajectory(s_start, s_mid,N_m, params)
%INIT_TRAJECTORY computes the initial trajectory for the iterative process
% arguments: 
% s_start - starting point, the point at which the quad starts.
% N_m     - the amount of trajectory points at the mth stage.
% s_mid   - the middle point between the communication user and the
% estimated sensing target (the estimated position obtained from the
% radar).
% V_str   - the speed of the quad. It is constant at throughout the
% initial trajectory.
% T_f     - flight time between points.
% output:
% S_init  - the initial trajectory.
% V_out   - the velocity throughout the trajectory.

V_str = params.sim.V_str;
T_f   = params.sim.T_f;

S_init = zeros(2,N_m);
V_init  = ones(2,N_m).*V_str.*(s_mid - s_start)./norm(s_mid - s_start);
%S_init(1,:) = linspace(s_start(1),s_mid(1), N_m);
%S_init(2,:) = linspace(s_start(2),s_mid(2), N_m);

for n = 1:N_m
         % here we multiply by T_f. In the paper they dont, check that point
         S_init(:,n) = s_start + V_str .*n .* (s_mid - s_start) ./ norm(s_mid - s_start);
end

