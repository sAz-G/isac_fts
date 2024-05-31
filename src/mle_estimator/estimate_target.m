%------------------------------------------------------------------------
% FUNCTION NAME: estimate_target
% AUTHOR: Sharif Azem (sAz-G on github), Markus Krantzik (mardank on github)
%
% DESCRIPTION: estimate_target provides the estimated position of the target given the
% hover points. It uses one of the implemented algorithms to caluclate the
% maximum likelihood estimator.
%
% INPUTS:
%   S_hover - Hover points
%   D_meas  - Distance vector
%   params  - predefined parameters
%   algo    - type of algorithm to use (currently only grid search)
%
% OUTPUTS:
%   pos - x,y coordinates
%
% USAGE:  pos = estimate_target(S_hover,D_meas,params, algo)
%
%------------------------------------------------------------------------

function pos = estimate_target(S_hover,D_meas,params, algo)

if strcmp(algo,'gridsearch')
    % create grid vectors 
    [x_t,y_t] = grid_vectors(1500,1500,1000,1000); % create vectors for grid search
    % get x,y hover positions
    s_hover_x = S_hover(1,:);
    s_hover_y = S_hover(2,:);
    
    % get the estimated target matrix index
    [x_t_idx,y_t_idx]  = get_min(D_meas,x_t,y_t,s_hover_x,s_hover_y,params);
    pos = [x_t(x_t_idx);y_t(y_t_idx)];        
elseif strcmp(algo,'random_gridsearch')
    % create grid vectors
    rpts = 10;
    rng(3);
    X_t = 1500*rand(rpts,1000);
    rng(4);
    Y_t = 1500*rand(rpts,1000);
    
    % get x,y hover points
    s_hover_x = S_hover(1,:);
    s_hover_y = S_hover(2,:);
    
    % calculation of argmin
    min_prev = 1e12;
    for z = 1:rpts % iteration over repetitions
        x_t = X_t(z,:);
        y_t = Y_t(z,:);
        [x_t_idx,y_t_idx,min_val]  = get_min(D_meas,x_t,y_t,s_hover_x,s_hover_y,params);
        if min_val < min_prev % if the current min value is less than the previous
            min_prev = min_val;
            pos = [x_t(x_t_idx);y_t(y_t_idx)];
        end    
    end

end

end

