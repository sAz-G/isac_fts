function plot_pdf(X_t, Y_t, X_jt, Y_jt, H_mat, D_mat_est, params)
%PLOT_MLE Summary of this function goes here
%   Detailed explanation goes here

D_xy = sqrt((X_t-X_jt).^2+(Y_t-Y_jt).^2 + H_mat.^2);
diff_square = (D_mat_est - D_xy).^2;

P              = params.sim.P;
G_p            = params.sim.G_p;
beta_0         = params.sim.beta_0;
a              = params.sim.a;
sigma_0        = params.sim.sigma_0;

factor = (a * sigma_0.^2)./(P*G_p*beta_0);
sigma_m_square = factor .* D_xy.^4;
exp_factor = 1./sqrt(2*pi .* sigma_m_square);

before_product = exp_factor .* exp(-1./(2*sigma_m_square) .* diff_square);

func = prod(before_product, 3);

figure
[~,idx] = max(func(:));
[idy,idx] = ind2sub(size(func),idx);

x_range = [0, params.sim.L_x];
y_range = [0, params.sim.L_y];

imagesc(x_range, y_range, func);
set(gca,'YDir','normal')

xlabel("x in m")
ylabel("y in m")

% hold on;
% plot(idx * params.sim.L_x/length(func(:,1)), idy * params.sim.L_y/length(func(1,:)), "*r");

title("PDF of the position of the target");

end

