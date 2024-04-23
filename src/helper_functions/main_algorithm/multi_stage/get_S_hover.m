%------------------------------------------------------------------------
% FUNCTION NAME: get_S_hover
% AUTHORS: Sharif Azem (sAz-G on GitHub) and Markus Krantzik (mardank on GitHub)
%
% DESCRIPTION: 
%   This function extracts hover points from a given trajectory S within
%   the specified range of stages [m_start, m_end]. Hover points represent
%   locations where the trajectory remains stationary before advancing to
%   the next stage.
%
% INPUTS:
%   S         - Trajectory at stage m, a 2-by-N matrix representing the x
%               and y coordinates of N points along the trajectory.
%   params    - Structure containing simulation and setup parameters,
%               including 'sim' which contains N_stg and mu.
%   m_start   - Starting stage number.
%   m_end     - Ending stage number.
%
% OUTPUTS:
%   hover_idxs - Index vector of the hover points.
%   S_hover    - Hover points up to the (m-1)th stage, a 2-by-M matrix
%                containing the x and y coordinates of M hover points.
%
% USAGE: 
%   [hover_idxs, S_hover] = get_S_hover(params, m_start, m_end, S)
%
%------------------------------------------------------------------------

function [hover_idxs, S_hover] = get_S_hover(params, m_start, m_end, S)

    % Extract simulation parameters
    N_stg = params.sim.N_stg;
    mu    = params.sim.mu;
    K_stg = floor(N_stg / mu);

    % Calculate the shift using modulo
    shft        = mod(N_stg, mu); 
    shft_vec    = shft * mod((m_start - 1 + mu):(m_end - 1 + mu), mu);
    shft_mat    = ones(K_stg, m_end - m_start + 1) .* shft_vec;
    idxs_shift  = shft_mat(:)'; 

    % Calculate hover indices
    hover_idxs_linear   = mu + (m_start - 1) * K_stg * mu : mu : (m_end * K_stg * mu);
    hover_idxs = idxs_shift + hover_idxs_linear;

    % Extract hover points
    if ~isempty(S)
        S_hover = S(:, hover_idxs);
    else
        S_hover = [];
    end

end
