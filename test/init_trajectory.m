function [S_out, V_out] = init_trajectory(s_start, N_m, x_mid, V_str, T_f)
%INIT_TRAJECTORY Summary of this function goes here
%   Detailed explanation goes here
S_out = zeros(2,N_m);
V_out = zeros(2,N_m);

for n = 1:N_m
        S_out(:,n) = s_start + V_str * n * (x_mid - s_start) / norm(x_mid - s_start);

        if n == 1
            V_out(:,n) = (S_out(:, n) - s_start) / T_f;
        else
            V_out(:,n) = (S_out(:, n) - S_out(:, n-1)) / T_f;
        end
end

