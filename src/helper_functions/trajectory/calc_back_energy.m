function E_back = calc_back_energy(S_start, s_base, params)
%CALC_BACK_ENERGY Summary of this function goes here
%   Detailed explanation goes here
N_back = ceil(norm(S_start-s_base) ./ params.sim.V_max / params.sim.T_f);
K_back = floor(N_back/params.sim.mu);

%V_back = [ones(1,N_back)*params.sim.V_max; zeros(1,N_back)];
S_back = [linspace(S_start(1), s_base(1), N_back+1); linspace(S_start(2), s_base(2), N_back+1)];
S_back = S_back(:, 2:end);

E_back = calc_real_energy(K_back, S_back, S_start, params);
end

