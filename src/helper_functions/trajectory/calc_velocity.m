function V = calc_velocity(S, S_s, params)
%   calc_velocity Summary of this function goes here
%   Detailed explanation goes here
V = diff([S_s, S], 1, 2)/params.sim.T_f;
end

