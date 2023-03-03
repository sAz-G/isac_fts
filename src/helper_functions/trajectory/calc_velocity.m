function V = calc_velocity(S, s_s, params)
%   calc_velocity Summary of this function goes here
%   Detailed explanation goes here
V = diff([s_s, S], 1, 2)./params.sim.T_f;
end

