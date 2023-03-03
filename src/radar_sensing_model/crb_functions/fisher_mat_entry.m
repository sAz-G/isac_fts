function fisher_entry = fisher_mat_entry(S_hov, s_target,entry,params,relative_dist_vec,factor_CRB)
%FISHER_MAT_ENTRY Summary of this function goes here
%   Detailed explanation goes here

if ~exist('relative_dist_vec', 'var')
        relative_dist_vec         = relative_distance(S_hov, s_target, params.sim.H);
end

if ~exist('factor_CRB', 'var')
        P              = params.sim.P;
        G_p            = params.sim.G_p;
        beta_0         = params.sim.beta_0;
        a              = params.sim.a;
        sigma_0        = params.sim.sigma_0;
        factor_CRB     = (P.*G_p.*beta_0)/(a*sigma_0.^2);
end

if strcmp(entry,'theta_a')    
    relative_pos_vec_x = (S_hov(1,:) - s_target(1))';
    
    left_sumand  = factor_CRB.*(relative_pos_vec_x.'*diag(1./relative_dist_vec.^6)*relative_pos_vec_x);
    right_sumand = 8.*(relative_pos_vec_x.'*diag(1./relative_dist_vec.^4)*relative_pos_vec_x);
   
    fisher_entry = left_sumand + right_sumand;

elseif strcmp(entry,'theta_b')
    fisher_entry =...
        fisher_mat_entry([S_hov(2,:);S_hov(1,:)], [s_target(2);s_target(1)],...
                         'theta_a',params,relative_dist_vec);
elseif strcmp(entry,'theta_c')
    relative_pos_x = (S_hov(1,:) - s_target(1)).';
    relative_pos_y = (S_hov(2,:) - s_target(2)).';

    fisher_entry = entry_c(relative_pos_x, relative_pos_y, relative_dist_vec,factor_CRB);
    
elseif strcmp(entry,'all')
else
end

function theta_c = entry_c(relative_pos_x, relative_pos_y, relative_dist_vec,factor_CRB)
    
    left  = factor_CRB.*(relative_pos_x.'*diag(1./relative_dist_vec.^6)*relative_pos_y);
    right = 8.*(relative_pos_x.'*diag(1./relative_dist_vec.^4)*relative_pos_y);
   
    theta_c = left + right;
end


end

