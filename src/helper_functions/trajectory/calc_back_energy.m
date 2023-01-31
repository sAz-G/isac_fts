function E_back = calc_back_energy(s_q, s_b, params)
%CALC_BACK_ENERGY Summary of this function goes here
%   Detailed explanation goes here
N_back = ceil(norm(s_q-s_b)./params.sim.V_max);
K_back = floor(N_back/params.sim.mu);
V_back = [ones(1,N_back)*params.sim.V_max; zeros(1,N_back)];

E_back = calc_real_energy(K_back,V_back,params);
end

