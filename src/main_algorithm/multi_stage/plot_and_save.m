clear
close all

folder_structure = ls;
num_elements = length(folder_structure(:,1));

for i = 3:num_elements
    if contains(folder_structure(i,:),'.mat')
        load(folder_structure(i,:))
        break
    end
end

plot_map(S_opt_mat, s_b, s_t, S_target_est_mat, s_c, params)

file_name = workspace_name + '_map';
saveas(gcf, file_name, 'epsc');
saveas(gcf, file_name, 'png');

figure
plot(res.J(~isnan(res.J(:,1)),~isnan(res.J(1,:)))');

title('Objective function over the iterations per stage')
ylabel('$\widetilde{\mathrm{CRB}}_\mathrm{Taylor}^m - (1 - \eta) \cdot \overline{R}^m_{\mathrm{Taylor}}$','interpreter','latex')
xlabel('iterations')

dim = size(res.J(~isnan(res.J(:,1)),~isnan(res.J(1,:))));

legend_label = strings(dim(1),1);
for i = 1:dim(1)
    legend_label(i) = 'm=' + string(i);
end

l = legend(legend_label);
l.Location = 'southeast';

grid on

file_name = workspace_name + '_oscillation';
saveas(gcf, file_name, 'epsc');
saveas(gcf, file_name, 'png');
