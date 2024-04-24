%------------------------------------------------------------------------
% FUNCTION NAME: get_min_gridsearch
% AUTHORS: Sharif Azem (sAz-G on GitHub), Markus Krantzik (mardank on GitHub)
%
% DESCRIPTION: 
%   We want to find the minimum (objective function (of)).
%
% INPUTS:
%   x_t, y_t: Parameters of the MLE. These are the possible target coordinates.
%   D_est: The vector of estimated distances.
%   x_jt, y_jt: Coordinates of the hovering positions.
%   params: Parameters for the calculation.
%
% OUTPUTS:
%   idx: Matrix id of the x position.
%   idy: Matrix id of the y position.
%   min_val: Minimum value.
%
% USAGE: 
%   [idx, idy, min_val] = get_min_gridsearch(D_est, x_t, y_t, x_jt, y_jt, params)
%
%------------------------------------------------------------------------

function [idx, idy, min_val] = get_min_gridsearch(D_est, x_t, y_t, x_jt, y_jt, params)

    x_dim = length(x_t);    % Dimension in x direction of the matrices after using meshgrid
    y_dim = length(y_t);    % Dimension in y direction of the matrices after using meshgrid
    xj_dim = length(x_jt);  % Dimension in x direction of the matrices after using meshgrid
    yj_dim = length(y_jt);  % Dimension in y direction of the matrices after using meshgrid

    % Create a grid
    [X_t, Y_t] = meshgrid(x_t, y_t);
    X_t = repmat(X_t, [1,1,xj_dim]);
    Y_t = repmat(Y_t, [1,1,yj_dim]);

    % Create distance grid and grid of hovering positions
    D_mat_est = ones(x_dim, y_dim, length(D_est));
    Y_jt = ones(x_dim, y_dim, yj_dim);
    X_jt = ones(x_dim, y_dim, xj_dim);
    for k = 1:length(D_est)
        D_mat_est(:,:,k) = D_est(k);
        Y_jt(:,:,k) = y_jt(k);
        X_jt(:,:,k) = x_jt(k);
    end

    % Calculate the MLE
    P_xy = (X_t - X_jt).^2 + (Y_t - Y_jt).^2 + params.sim.H.^2;
    log_Pxy = log(P_xy);
    numerator = D_mat_est - sqrt(P_xy);
    denominator = P_xy;
    fraction = (numerator ./ denominator).^2;

    % Constants
    P = params.sim.P;
    G_p = params.sim.G_p;
    beta_0 = params.sim.beta_0;
    a = params.sim.a;
    sigma_0 = params.sim.sigma_0;

    % Final expression of the MLE
    factor = (P * G_p * beta_0) / (a * sigma_0.^2);
    func = sum(log_Pxy + factor * fraction, 3);

    % Calculate min 
    [min_val, idx] = min(func(:));
    [idy, idx] = ind2sub(size(func), idx); 

end
