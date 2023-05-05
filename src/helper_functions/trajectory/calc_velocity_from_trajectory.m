function V_out = calc_velocity_from_trajectory(S, S_s, params)
%CALC_VELOCITY_FROM_TRAJECTORY Calulate the velocity from a given
%trajectory of the drone
%   S == trajectory of the flight
%   S_s == start point of the trajectory
%   params == struct with simulation parameters

V_out = diff([S_s, S], 1, 2)/params.sim.T_f;
end

