function S_out = init_trajectory(s_start, N_m, x_mid, V_str)
%INIT_TRAJECTORY Summary of this function goes here
%   Detailed explanation goes here
S_out = zeros(2,N_m);

for i = 1:N_m
        S_out(:,i) = s_start + V_str*i*(x_mid - s_start) / norm(x_mid - s_start);
end

