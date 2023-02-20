function grad = fisher_entry_gradient(S_hov,s_t,type, direction,params,relative_dist_vec,factor_CRB)
%FISHER_ENTRY_TAYLOR Summary of this function goes here
%   Detailed explanation goes here

if ~exist('relative_dist_vec', 'var')
        relative_dist_vec         = relative_distance(S_hov, s_t, params.sim.H);
end

if ~exist('factor_CRB', 'var')
        P              = params.sim.P;
        G_p            = params.sim.G_p;
        beta_0         = params.sim.beta_0;
        a              = params.sim.a;
        sigma_0        = params.sim.sigma_0;
        factor_CRB     = (P.*G_p.*beta_0)/(a*sigma_0.^2);
end


if strcmp(type, 'theta_a')
    rel_pos_x = S_hov(1,:) - s_t(1);    
    rel_pos_y = S_hov(2,:) - s_t(2);   

    if     strcmp(direction, 'x')
        grad = factor_CRB.*(2.*rel_pos_x.*relative_dist_vec.^2-6.*rel_pos_x.^3)./relative_dist_vec.^8+ ...
               (16.*rel_pos_x.*relative_dist_vec.^2-24.*rel_pos_x.^3)./relative_dist_vec.^6;
    elseif strcmp(direction, 'y')
        grad = rel_pos_y.*rel_pos_x.^2.*(24./relative_dist_vec.^6 ...
                                        -6.*factor_CRB./relative_dist_vec.^8);
    end
elseif  strcmp(type, 'theta_b')

    if  strcmp(direction, 'x')
       grad =  fisher_entry_gradient([S_hov(2,:);S_hov(1,:)], [s_t(2);s_t(1)],'theta_a', 'y',...
                               params,relative_dist_vec,factor_CRB);
    elseif strcmp(direction, 'y')
       grad = fisher_entry_gradient([S_hov(2,:);S_hov(1,:)], [s_t(2);s_t(1)],'theta_a', 'x',...
                               params,relative_dist_vec,factor_CRB);
    end
    
elseif  strcmp(type, 'theta_c')
    grad = c_entry_grad(S_hov, s_t, relative_dist_vec,direction,factor_CRB);
else

end


function grad_c = c_entry_grad(S_hov, s_t, relative_dist_vec,direction,factor_CRB)
    if     strcmp(direction, 'x')
            re_pos_x = S_hov(1,:) - s_t(1);    
            re_pos_y = S_hov(2,:) - s_t(2);   
    elseif strcmp(direction, 'y')
            re_pos_y = S_hov(1,:) - s_t(1);    
            re_pos_x = S_hov(2,:) - s_t(2);   
    end
    grad_c = re_pos_y./(relative_dist_vec.^6).*(8.*relative_dist_vec.^2 - 24.*re_pos_x.^2 ...
                                                 +factor_CRB./relative_dist_vec.^2.*(relative_dist_vec.^2-6.*re_pos_x.^2));
end

end

