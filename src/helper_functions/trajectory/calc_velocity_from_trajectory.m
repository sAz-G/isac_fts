function V_out = calc_velocity_from_trajectory(S, S_s, params)
%CALC_VELOCITY_FROM_TRAJECTORY Summary of this function goes here
%   Detailed explanation goes here
V_out = diff([S_s, S], 1, 2)/params.sim.T_f;
end

