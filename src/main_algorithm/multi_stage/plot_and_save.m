close all

plot_map(S_opt_mat, s_b, s_t, S_target_est_mat, s_c, params)

file_name = workspace_name + '_map';
saveas(gcf, fullfile(out_path, file_name), 'epsc');
saveas(gcf, fullfile(out_path, file_name), 'png');

figure
plot(res.J(~isnan(res.J(:,1)),~isnan(res.J(1,:)))');

title('Objective function over the iterations per stage')
ylabel('$\widetilde{\mathrm{CRB}}_\mathrm{Taylor}^m - (1 - \eta) \cdot \overline{R}^m_{\mathrm{Taylor}}$','interpreter','latex')
xlabel('iterations')

file_name = workspace_name + '_oscillation';
saveas(gcf, fullfile(out_path, file_name), 'epsc');
saveas(gcf, fullfile(out_path, file_name), 'png');

l = legend('m=1','m=2','m=3','m=4','m=5');
l.Location = 'southeast';

grid on