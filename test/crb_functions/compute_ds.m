function ds_out = compute_ds(S_hover_in, S_target_est_in, H_in, K_tot_in)
%COMPUTE_GRADIENT Summary of this function goes here
%   Detailed explanation goes here
    x_target_est_diff = (S_hover_in(1,:) - S_target_est_in(1));
    y_target_est_diff = (S_hover_in(2,:) - S_target_est_in(2));
    
    % d_s_km = sqrt(H^2 + x_target_est_diff.^2 + y_target_est_diff.^2);
    % d_s_km = norms([x_target_est_diff; y_target_est_diff; H_in*ones(1, K_tot_in)], 2, 1);
    ds_out = norms([x_target_est_diff; y_target_est_diff; H_in*ones(1, K_tot_in)], 2, 1);
end

