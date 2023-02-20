function fisher_entry = fisher_mat_entry(S_hov, s_target,entry,params,relative_dist_vec)
%FISHER_MAT_ENTRY Summary of this function goes here
%   Detailed explanation goes here

if ~exist('relative_dist_vec', 'var')
        relative_dist_vec         = relative_distance(S_hov, s_target, params.sim.H);
end

if strcmp(entry,'theta_a')    
    relative_pos_vec_x = (S_hov(1,:) - s_target(1))';
    
    P              = params.sim.P;
    G_p            = params.sim.G_p;
    beta_0         = params.sim.beta_0;
    a              = params.sim.a;
    sigma_0        = params.sim.sigma_0;
    factor_CRB     = (P.*G_p.*beta_0)/(a*sigma_0.^2);
    
    left_sumand  = factor_CRB.*(relative_pos_vec_x.'*diag(1./relative_dist_vec.^6)*relative_pos_vec_x);
    right_sumand = 8.*(relative_pos_vec_x.'*diag(1./relative_dist_vec.^4)*relative_pos_vec_x);
   
    fisher_entry = left_sumand + right_sumand;

elseif strcmp(entry,'theta_b')
    fisher_entry =...
        fisher_mat_entry([S_hov(2,:);S_hov(1,:)], [s_target(2);s_target(1)],...
                         'theta_a',params,relative_dist_vec);
elseif strcmp(entry,'theta_c')
    
    x_entries = sqrt((S_hov(2,:)-s_target(2)).*(S_hov(1,:)-s_target(1)));
    y_entries = sqrt((S_hov(1,:)-s_target(1)).^2 + (S_hov(2,:)-s_target(2)).^2 - x_entries.^2);
    
    fisher_entry =...
        fisher_mat_entry([x_entries;y_entries], [0,0],...
                        'theta_a',params,relative_dist_vec);

elseif strcmp(entry,'all')
else
end

end

