function relative_distance = user_quad_distance(S,s_c,H)
%USER_QUAD_DISTANCE Summary of this function goes here
%   Detailed explanation goes here
    relative_pos_x = (S(1,:) - s_c(1));
    relative_pos_y = (S(2,:) - s_c(2));
    
    N = size(S,2);
    relative_distance = norms([relative_pos_x; relative_pos_y; H*ones(1, N)], 2, 1);
end

