%------------------------------------------------------------------------
% FUNCTION NAME: get_S_hover
% AUTHORS: Sharif Azem (sAz-G on GitHub), Markus Krantzik (mardank on GitHub)
%
% DESCRIPTION: Extract the hover points out of a given trajectory S given
%              range of stages [m_start, m_end].
%
% INPUTS:
%   S         - The trajectory S at the stage m
%   params    - Simulation and setup parameters
%   m_start   - Starting stage number  
%   m_end     - Ending stage number 
%
% OUTPUTS:
%   S_hover   - Hover points (already visited) up to the (m-1)th stage
%   hover_idx - Index vector of the hover points 
%
% USAGE: 
%   [hover_idxs, S_hover] = get_S_hover(params, m_start, m_end, S)
%
%------------------------------------------------------------------------

function [hover_idxs, S_hover] = get_S_hover(params, m_start, m_end, S)

    N_stg = params.sim.N_stg;
    mu    = params.sim.mu;
    K_stg = floor(N_stg ./ mu);

    % Calculate the shift using modulo
    shft        = mod(N_stg, mu); 
    % Shift for all stages 
    shft_vec    = shft * mod((m_start - 1 + mu):(m_end - 1 + mu), mu);
    shft_mat    = ones(K_stg, m_end - m_start + 1) .* shft_vec;
    idxs_shift  = shft_mat(:)'; 
    % Linear hovering indices
    hover_idxs_linear = mu + (m_start - 1) * K_stg * mu : mu : (m_end * K_stg * mu);
    % Multi-stage hovering indices
    hover_idxs = idxs_shift + hover_idxs_linear;

    % Extract the hover points after calculating the index
    if ~isempty(S)
        S_hover = S(:, hover_idxs);
    else
        S_hover = [];
    end

end
