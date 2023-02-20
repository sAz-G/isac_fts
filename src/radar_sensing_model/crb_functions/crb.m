function crb_val = crb(S, s_t, params)
%CRB Summary of this function goes here
%   Detailed explanation goes here
fisher_a = fisher_mat_entry(S,s_t,'theta_a',params);
fisher_b = fisher_mat_entry(S,s_t,'theta_b',params);
fisher_c = fisher_mat_entry(S,s_t,'theta_c',params);

determinant = fisher_a.*fisher_b - fisher_c.^2;

crb_val = (fisher_a + fisher_b)./determinant;

end

