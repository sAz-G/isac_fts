%------------------------------------------------------------------------
% FUNCTION NAME: get_S_hover
% AUTHOR: Sharif Azem     (TU-Darmstadt department 18, sAz-G on github)
%         Markus Krantzik (TU-Darmstadt department 18, mardank on github)
%
% DESCRIPTION: extract the hover points out of a given trajectory S given
%              range of stages [m_start, m_end].
%
% INPUTS:
% S         - the trajectory S at the stage m
% params    - simulation and setup parameters
% m_start   - starting stage number  
% m_end     - ending stage number 
%
% OUTPUTS:
% S_hover   - hover points (already visited) up to the (m-1)th stage
% hover_idx - index vector of the hover points 
%
% USAGE: [hover_idxs, S_hover] = get_S_hover( params, m_start,m_end,S)
%
%------------------------------------------------------------------------

function [hover_idxs, S_hover] = get_S_hover(params, m_start,m_end,S)

N_stg = params.sim.N_stg;
mu    = params.sim.mu;
K_stg = floor(N_stg./mu);

% calc the shift using modulo
shft                = mod(N_stg, mu); 
% shift for all stages 
shft_vec            = shft*mod((m_start-1+mu):(m_end-1 + mu), mu);
shft_mat            = ones(K_stg, m_end-m_start+1).*shft_vec;
idxs_shift          = shft_mat(:)'; 
%linear hovering indices
hover_idxs_linear   = mu+(m_start-1)*K_stg*mu:mu:(m_end*K_stg*mu);
% multi stage hovering indices
hover_idxs          = idxs_shift + hover_idxs_linear;

% extract the hover points after calculating the index
if ~isempty(S)
    S_hover             = S(:,hover_idxs);
else
    S_hover = [];
end

end

