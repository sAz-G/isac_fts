function [hover_idxs, S_hover] = get_S_hover( params, m_start,m_end,S)
%GET_S_HOVER extract the hover points out of a given trajectory S given
%   range of stages [m_start, m_end].
% The hover points of the trajectory S when using the multistage algorithm differ
% from the hover points when calculation the optimal trajectory one time
% (S_hover_multi_stage(1:end) != S_hover_one_stage, with
% size(S_hover_multi_stage) == size(S_hover_one_stage).
%
% Input: 
% S         - the trajectory S at the stage m
% params    - simulation and setup parameters
% m_start   - starting stage number  
% m_end     - ending stage number 
%
% Output:
% S_hover   - hover points (already visited) up to the (m-1)th stage
% hover_idx - index vector of the hover points 


N_stg = params.sim.N_stg;
mu    = params.sim.mu;
K_stg = floor(N_stg./mu);

shft                = mod(N_stg, mu); 
shft_vec            = shft*mod((m_start-1+mu):(m_end-1 + mu), mu);
shft_mat            = ones(K_stg, m_end-m_start+1).*shft_vec;
idxs_shift          = shft_mat(:)'; 
hover_idxs_linear   = [mu+(m_start-1)*K_stg*mu:mu:(m_end*K_stg*mu)];
hover_idxs          = idxs_shift + hover_idxs_linear;

if ~isempty(S)
    S_hover             = S(:,hover_idxs);
else
    S_hover = [];
end

end

