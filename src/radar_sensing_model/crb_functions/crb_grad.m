function gradient_vector = crb_grad(S,s_t,params,dimension)
%CRB_TAYLOR Summary of this function goes here
%   Detailed explanation goes here

fisher_a = fisher_mat_entry(S,s_t,'theta_a',params);
fisher_b = fisher_mat_entry(S,s_t,'theta_b',params);
fisher_c = fisher_mat_entry(S,s_t,'theta_c',params);

fisher_determinant  = fisher_a.*fisher_b - fisher_c.^2;

K_stg       = floor(params.sim.N_stg./params.sim.mu);
S_grad    = S(:,end-K_stg+1:end);

fisher_a_grad_dim = fisher_entry_gradient(S_grad,s_t,'theta_a', dimension,params);
fisher_b_grad_dim = fisher_entry_gradient(S_grad,s_t,'theta_b', dimension,params);
fisher_c_grad_dim = fisher_entry_gradient(S_grad,s_t,'theta_c', dimension,params);

fisher_det_grad_dim = (fisher_a_grad_dim.*fisher_b + fisher_b_grad_dim.*fisher_a - 2.*fisher_c_grad_dim.*fisher_c);

denominator = fisher_b_grad_dim.*fisher_determinant -...
              fisher_b.*fisher_det_grad_dim;

denominator = denominator + fisher_a_grad_dim.*fisher_determinant -...
                            fisher_a.*fisher_det_grad_dim;

gradient_vector = denominator./fisher_determinant.^2;

end

