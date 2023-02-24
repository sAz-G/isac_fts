function fig_conv = plot_convergence(debug_S, debug_V, debug_delta, debug_xi, S_s, m, params)
%PLOT_CONVERGENCE Summary of this function goes here
%   Detailed explanation goes here

fig_conv = figure();
iter = params.sim.iter;
T_f = params.sim.T_f;

title("Convergence of the optimization variable at Stage " + m)

subplot(4, 2, 1);
title("UAV trajectory over the iterations")
xlabel("x in m")
ylabel("y in m")
hold on;
grid on;

subplot(4, 2, 2);
title("Mean velocity over the iterations")
xlabel("iterations")
ylabel("Mean velocity")
hold on;
grid on;

subplot(4, 2, 3);
title("Mean of delta over the iterations")
xlabel("iterations")
ylabel("Mean delta")
hold on;
grid on;

subplot(4, 2, 4);
title("Mean of xi over the iterations")
xlabel("iterations")
hold on;
grid on;
xlabel("iterations")
ylabel("Mean xi")

subplot(4, 2, [5, 6]);
title("Difference between delta^2 and xi")
xlabel("iterations")
hold on;
grid on;
xlabel("iterations")
ylabel("Difference between delta^2 and xis")

subplot(4, 2, [7, 8]);
title("Velocity in S and V")
xlabel("iterations")
hold on;
grid on;
xlabel("iterations")
ylabel("S an V")

for i = 1:iter
    subplot(4, 2, 1);
    if i ~= iter
        plot(debug_S(1,:,i), debug_S(2,:,i), "k");
    else
        plot(debug_S(1,:,i), debug_S(2,:,i), "r");
    end
    
    subplot(4, 2, 2);
    plot(i, mean(norms(debug_V(:,:,i), 2, 1)), "ob");

    subplot(4, 2, 3);
    plot(i, mean(debug_delta(:,:,i)), "ob");

    subplot(4, 2, 4);
    plot(i, mean(debug_xi(:,:,i)), "ob");

    subplot(4, 2, [5, 6]);
    hold on;
    plot(i, mean(abs(debug_xi(:,:,i) - debug_delta(:,:,i).^2)), "ob")
    
    subplot(4, 2, [7, 8]);
    hold on;
    plot(i, mean(norms([diff(debug_S(1,:,1)); diff(debug_S(2,:,1))]./T_f, 2, 1)), "ob");
    plot(i, mean(norms(debug_V(:,2:end,1), 2, 1)), "or");
    
    
end
subplot(4, 2, [7, 8]);
hold on;
legend("V from S", "V from V")

title("Parameter Convergence at Stage " + m)

end

