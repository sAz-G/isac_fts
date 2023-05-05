function dc = compute_dc(S_target_in, S_c, H_in, N_in)
%COMPUTE_DC Compute a vector of distances to communication user over a trajectory 
%   S_target_in == Position of the drone
%   S_c == position of the communication user
%   H_in == flying height of the drone
%   N_in == number if set_points
%   params == struct with simulation parameters

    x_c_diff = (S_target_in(1,:) - S_c(1));
    y_c_diff = (S_target_in(2,:) - S_c(2));
    
    % distance vector
    dc = norms([x_c_diff; y_c_diff; H_in*ones(1, N_in)], 2, 1);
end

