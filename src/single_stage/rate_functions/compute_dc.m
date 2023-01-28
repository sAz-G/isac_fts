function dc = compute_dc(S_target_in, S_c, H_in, N_in)
%COMPUTE_DC Summary of this function goes here
%   Detailed explanation goes here
    x_c_diff = (S_target_in(1,:) - S_c(1));
    y_c_diff = (S_target_in(2,:) - S_c(2));
    
    dc = norms([x_c_diff; y_c_diff; H_in*ones(1, N_in)], 2, 1);
end

