%------------------------------------------------------------------------
% FUNCTION NAME: crb_grad
% AUTHORS: Sharif Azem (sAz-G on GitHub), Markus Krantzik (mardank on GitHub)
%
% DESCRIPTION: 
%   Calculate the gradient of the CRB function.
%
% INPUTS:
%   S - Hover points of the drone.
%   s_t - Position of the sensing target.
%   params - Predefined parameters.
%   dimension - The dimension, along which the gradient is calculated (x or y).
%
% OUTPUTS:
%   gradient_vector - A vector of the gradient.
%
% USAGE: 
%   gradient_vector = crb_grad(S, s_t, params, dimension);
%
%------------------------------------------------------------------------

function gradient_vector = crb_grad(S, s_t, params, dimension)

    % Fisher entries
    fisher_a = fisher_mat_entry(S, s_t, 'theta_a', params);
    fisher_b = fisher_mat_entry(S, s_t, 'theta_b', params);
    fisher_c = fisher_mat_entry(S, s_t, 'theta_c', params);

    % Determinant 
    fisher_determinant  = fisher_a .* fisher_b - fisher_c.^2;

    % Points, for which the gradient is to be calculated
    K_stg = params.sim.K_stg;   
    S_grad = S(:, end - K_stg + 1:end); 

    % Gradient of the Fisher entry
    fisher_a_grad_dim = fisher_entry_gradient(S_grad, s_t, 'theta_a', dimension, params);
    fisher_b_grad_dim = fisher_entry_gradient(S_grad, s_t, 'theta_b', dimension, params);
    fisher_c_grad_dim = fisher_entry_gradient(S_grad, s_t, 'theta_c', dimension, params);

    % Gradient of the Fisher determinant
    fisher_det_grad_dim = (fisher_a_grad_dim .* fisher_b + fisher_b_grad_dim .* fisher_a - 2 .* fisher_c_grad_dim .* fisher_c);

    % Calculate the final result
    denominator = fisher_b_grad_dim .* fisher_determinant -...
                  fisher_b .* fisher_det_grad_dim;

    denominator = denominator + fisher_a_grad_dim .* fisher_determinant -...
                                  fisher_a .* fisher_det_grad_dim;

    gradient_vector = denominator ./ fisher_determinant.^2;

end
