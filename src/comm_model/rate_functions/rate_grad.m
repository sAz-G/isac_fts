function grad_rate_dim = rate_grad(S,s_c,params,dim)
%RATE_GRADIENT_DIM Summary of this function goes here
%   Detailed explanation goes here

P              = params.sim.P;
sigma_0        = params.sim.sigma_0;
alpha_0        = params.sim.alpha_0;
B              = params.sim.B;
N              = size(S,2);
N_stg          = params.sim.N_stg;

rel_dist       = user_quad_distance(S(:,end-N_stg+1:end), s_c, params.sim.H);

if strcmp(dim,'x')
    rel_pos_x     = (S(1,end-N_stg+1:end) - s_c(1));
    grad_distance = rel_pos_x./rel_dist;
    inner_grad    = (P.*alpha_0)./(sigma_0.^2).*(-2)./(rel_dist.^3).*grad_distance;
    grad_rate_dim = 1./(1 + (P * alpha_0)./(sigma_0.^2.*rel_dist.^2)).*inner_grad.*B./(N * log(2));
elseif strcmp(dim,'y')
    rel_pos_y     = (S(2,end-N_stg+1:end) - s_c(2));
    grad_distance = rel_pos_y./rel_dist;
    inner_grad    = (P.*alpha_0)./(sigma_0.^2).*(-2)./(rel_dist.^3).*grad_distance;
    grad_rate_dim = 1./(1 + (P * alpha_0)./(sigma_0.^2*rel_dist.^2)).*inner_grad.*B./(N * log(2)) ;
end


end

