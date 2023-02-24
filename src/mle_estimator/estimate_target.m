function pos = estimate_target(S_hover,D_meas,params, algo)
% estimate_target provides the estimated position of the target given the
% hover points. It uses one of the implemented algorithm to caluclate the
% maximum likelihood estimator.
% 
% input: 
% S_hover - hover points. Points at which the quad measured the distance to
%           the target.
% D_meas  - distance measurement vector
% params  - predefined constant parameters 
% type    - type of output. (index, positions, both)
% algo    - which algorithm to use for the calculation of the maximum
%           likelihood estimator

if strcmp(algo,'gridsearch')
    [x_t,y_t] = grid_vectors(1500,1500,1000,1000); % create vectors for grid search
    s_hover_x = S_hover(1,:);
    s_hover_y = S_hover(2,:);
    
    % get the estimated target matrix index
    [x_t_idx,y_t_idx]  = get_min(D_meas,x_t,y_t,s_hover_x,s_hover_y,params);
    pos = [x_t(x_t_idx);y_t(y_t_idx)];        
elseif strcmp(algo,'random_gridsearch')
    % not implemmented yet
elseif strcmp(algo,'grad_dec')
    % not implemmented yet
elseif strcmp(algo,'cvx_app')
    % not implemmented yet
else
end

end
