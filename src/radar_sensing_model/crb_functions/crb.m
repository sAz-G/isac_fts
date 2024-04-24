%------------------------------------------------------------------------
% FUNCTION NAME: crb
% AUTHORS: Sharif Azem (sAz-G on GitHub), Markus Krantzik (mardank on GitHub)
%
% DESCRIPTION: 
%   Calculate the CRB of given points and target position.
%
% INPUTS:
%   S - Hover points of the drone.
%   s_t - Position of the sensing target.
%   params - Predefined parameters.
%
% OUTPUTS:
%   crb_val - CRB value.
%
% USAGE: 
%   crb_val = crb(S, s_t, params);
%
%------------------------------------------------------------------------

function crb_val = crb(S, s_t, params)
    fisher_a = fisher_mat_entry(S, s_t, 'theta_a', params);
    fisher_b = fisher_mat_entry(S, s_t, 'theta_b', params);
    fisher_c = fisher_mat_entry(S, s_t, 'theta_c', params);

    determinant = fisher_a .* fisher_b - fisher_c.^2;

    crb_val = (fisher_a + fisher_b) ./ determinant;
end
