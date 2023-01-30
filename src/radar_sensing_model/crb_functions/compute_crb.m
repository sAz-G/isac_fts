function CRB_out= compute_crb(S_traj_init_in, S_target_est_in, H_in, K_tot_in, params)
%COMPUTE_CRB Summary of this function goes here
%   Detailed explanation goes here
    
    [theta_a, theta_b, theta_c] = compute_thetas(S_traj_init_in, S_target_est_in, H_in, K_tot_in, params);

    CRB_out = (theta_b + theta_a)/(theta_a * theta_b - theta_c^2); % vectorized
end

