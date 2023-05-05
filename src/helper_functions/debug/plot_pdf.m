function plot_pdf(X_t, Y_t, X_jt, Y_jt, H_mat, D_mat_est, params)
%PLOT_MLE Plot the propability distrubution function

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
func_exp = prod(before_product, 3);

func_log = sum(log(exp_factor), 3) + sum(-1./(2*sigma_m_square) .* diff_square, 3);

figure
[~,idx] = max(func_exp(:));
[idy,idx] = ind2sub(size(func_exp),idx);

x_range = [0, params.sim.L_x];
y_range = [0, params.sim.L_y];

imagesc(x_range, y_range, func_exp);
set(gca,'YDir','normal')

xlabel("x in m")
ylabel("y in m")

hold on;
plot(idx * params.sim.L_x/length(func_exp(:,1)), idy * params.sim.L_y/length(func_exp(1,:)), "*r");

title("PDF of the position of the target");
% 
% test_idx = idy;
% test_idy = idx;
% 
% x_t = X_jt(test_idx,test_idy,:);
% x_t = x_t(:);
% y_t = Y_jt(test_idx,test_idy,:);
% y_t = y_t(:);
% 
% x_in = X_t(test_idx,test_idy,1);
% x_in = x_in(1);
% y_in = Y_t(test_idx,test_idy,1);
% y_in = y_in(1);
% 
% d_est = D_mat_est(test_idx,test_idy,:);
% d_est = d_est(:);
% 
% d_func = @(x, y) sqrt((x*ones(params.sim.K_stg, 1) - x_t).^2+(y*ones(params.sim.K_stg, 1) - y_t).^2 + (params.sim.H * ones(params.sim.K_stg, 1)).^2);
% 
% exp_factor_func = @(x, y) 1./sqrt(2*pi .* factor .* d_func(x, y).^4);
% pdf_func = @(x) -(prod(exp_factor_func(x(1), x(2)) .* vpa(exp(-1./(2*(factor .* d_func(x(1), x(2)).^4)) .* (d_est - d_func(x(1), x(2))).^2), 20)));
% 
% x_test = [x_in, y_in];
% pdf_func(x_test)
% 
% func_exp(test_idx, test_idy)
% 
% x0 = [750, 100]
% fminsearch(pdf_func, x0)

end

