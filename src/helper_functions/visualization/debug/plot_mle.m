function plot_mle(func, params)
%PLOT_MLE Summary of this function goes here
%   Detailed explanation goes here

figure
[~,idx] = min(func(:));
[idy,idx] = ind2sub(size(func),idx);

x_range = [0, params.sim.L_x];
y_range = [0, params.sim.L_y];
func_range = [min(func(:)), min(func(:)) * 1e2];

imagesc(x_range, y_range, func, func_range);
set(gca,'YDir','normal')

xlabel("x in m")
ylabel("y in m")

hold on;
plot(idx * params.sim.L_x/length(func(:,1)), idy * params.sim.L_y/length(func(1,:)), "*r");

title("MLE (minimum) of the position of the target");

end

