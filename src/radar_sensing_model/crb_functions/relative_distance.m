function ds_vec = relative_distance(S_hover, s_target, H)
%COMPUTE_GRADIENT Summary of this function goes here
%   Detailed explanation goes here
    relative_pos_x = S_hover(1,:) - s_target(1);
    relative_pos_y = S_hover(2,:) - s_target(2);
    
    ds_vec = norms([relative_pos_x; relative_pos_y; H*ones(size(relative_pos_x))], 2, 1);
end

