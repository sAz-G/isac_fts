%------------------------------------------------------------------------
% FUNCTION NAME: fisher_entry_gradient
% AUTHORS: Sharif Azem (sAz-G on GitHub), Markus Krantzik (mardank on GitHub)
%
% DESCRIPTION: 
%   Calculate the gradient of an entry of the Fisher information matrix.
%
% INPUTS:
%   S_hov - Hover points of the drone.
%   s_t - Position of the sensing target.
%   type - Type of entry: 'theta_a', 'theta_b', or 'theta_c'.
%   direction - Direction of the gradient: 'x' or 'y'.
%   params - Predefined parameters.
%   relative_dist_vec - Vector of the relative distances.
%   factor_CRB - A constant calculated using the sensing parameters.
%
% OUTPUTS:
%   grad - A vector of the gradient.
%
% USAGE: 
%   grad = fisher_entry_gradient(S_hov, s_t, type, direction, params, relative_dist_vec, factor_CRB);
%
%------------------------------------------------------------------------

function grad = fisher_entry_gradient(S_hov, s_t, type, direction, params, relative_dist_vec, factor_CRB)

    % Check if relative_dist_vec is inserted as input. Calculate if not.
    if ~exist('relative_dist_vec', 'var')
        relative_dist_vec = relative_distance(S_hov, s_t, params.sim.H);
    end

    % Check if factor_CRB is inserted. Calculate if not.
    if ~exist('factor_CRB', 'var')
        P = params.sim.P;
        G_p = params.sim.G_p;
        beta_0 = params.sim.beta_0;
        a = params.sim.a;
        sigma_0 = params.sim.sigma_0;
        factor_CRB = (P * G_p * beta_0) / (a * sigma_0.^2);
    end

    % Calculate the gradient
    if strcmp(type, 'theta_a')
        rel_pos_x = S_hov(1, :) - s_t(1);    
        rel_pos_y = S_hov(2, :) - s_t(2);   

        if strcmp(direction, 'x')
            grad = factor_CRB .* (2 .* rel_pos_x .* relative_dist_vec.^2 - 6 .* rel_pos_x.^3) ./ relative_dist_vec.^8 + ...
                   (16 .* rel_pos_x .* relative_dist_vec.^2 - 32 .* rel_pos_x.^3) ./ relative_dist_vec.^6;
        elseif strcmp(direction, 'y')
            grad = -rel_pos_y .* rel_pos_x.^2 .* (32 ./ relative_dist_vec.^6 ...
                                                 + 6 .* factor_CRB ./ relative_dist_vec.^8);
        end
    elseif strcmp(type, 'theta_b')
        if strcmp(direction, 'x')
            grad = fisher_entry_gradient([S_hov(2, :); S_hov(1, :)], [s_t(2); s_t(1)], 'theta_a', 'y', ...
                                          params, relative_dist_vec, factor_CRB);
        elseif strcmp(direction, 'y')
            grad = fisher_entry_gradient([S_hov(2, :); S_hov(1, :)], [s_t(2); s_t(1)], 'theta_a', 'x', ...
                                          params, relative_dist_vec, factor_CRB);
        end
    elseif strcmp(type, 'theta_c')
        grad = c_entry_grad(S_hov, s_t, relative_dist_vec, direction, factor_CRB);
    else
        % Handle other cases if needed
    end

    % A function to calculate the gradient of theta_c
    function grad_c = c_entry_grad(S_hov, s_t, relative_dist_vec, direction, factor_CRB)
        if strcmp(direction, 'x')
            re_pos_x = S_hov(1, :) - s_t(1);    
            re_pos_y = S_hov(2, :) - s_t(2);   
        elseif strcmp(direction, 'y')
            re_pos_y = S_hov(1, :) - s_t(1);    
            re_pos_x = S_hov(2, :) - s_t(2);   
        end
        grad_c = re_pos_y .* (factor_CRB ./ (relative_dist_vec.^6) + 8 ./ (relative_dist_vec.^4));
        grad_c = grad_c - re_pos_x.^2 .* re_pos_y .* (6 * factor_CRB ./ (relative_dist_vec.^8) + 32 ./ (relative_dist_vec.^6));
    end

end
