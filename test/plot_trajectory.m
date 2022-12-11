function plot_trajectory(S, X_s, X_t, X_c)
%PLOT_TRAJECTORY Summary of this function goes here
%   Detailed explanation goes here

figure
plot(S(1,:), S(2,:), "o-")
hold on
grid on

plot(X_s(1), X_s(2), "*")
plot(X_t(1), X_t(2), "*")
plot(X_c(1), X_c(2), "*")
end

