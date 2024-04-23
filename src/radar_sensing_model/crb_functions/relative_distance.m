%------------------------------------------------------------------------
% FUNCTION NAME: sigma_k
% AUTHOR: Sharif Azem     (TU-Darmstadt department 18, sAz-G on github)
%         Markus Krantzik (TU-Darmstadt department 18, mardank on github)
%
% DESCRIPTION: Calculate the relative distance to the sensing target
%
% INPUTS:
%   S_hover       - Hover points
%   s_target      - Target positions
%   H             - Height of the drone
%
% OUTPUTS:
%   ds_vec - vector of relative distances
%
% USAGE:   ds_vec = relative_distance(S_hover, s_target, H)
%
%------------------------------------------------------------------------

function ds_vec = relative_distance(S_hover, s_target, H)
    relative_pos_x = S_hover(1,:) - s_target(1);
    relative_pos_y = S_hover(2,:) - s_target(2);
    
    ds_vec = norms([relative_pos_x; relative_pos_y; H*ones(size(relative_pos_x))], 2, 1);
end

