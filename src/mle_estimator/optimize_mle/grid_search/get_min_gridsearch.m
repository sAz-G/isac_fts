%------------------------------------------------------------------------
% FUNCTION NAME: get_min
% AUTHOR: Sharif Azem
%         Markus Krantzik
%
% DESCRIPTION: %  We want to find the
% minimum(objective function (of)).
%
% INPUTS:
%   x_t,y_t: Parameters of the mle. These are the
%            possible target coordinates.
%   D_est:   The vector of estimated distances .
%
% OUTPUTS:
%   idx - matrix id of the x position 
%   idy - matrix id of the y position
%   min_val - minimum value 
%
% USAGE:  [idx, idy,min_val] = get_min_gridsearch(D_est,x_t,y_t,x_jt,y_jt,params)
%
%------------------------------------------------------------------------

function [idx, idy,min_val] = get_min_gridsearch(D_est,x_t,y_t,x_jt,y_jt,params)

x_dim = length(x_t);    % dimension in x direction of the matrices after using meshgrid
y_dim = length(y_t);    % dimension in y direction of the matrices after using meshgrid
xj_dim = length(x_jt);  % dimension in x direction of the matrices after using meshgrid
yj_dim = length(y_jt);  % dimension in y direction of the matrices after using meshgrid

% create a grid
[X_t, Y_t] = meshgrid(x_t, y_t);
X_t        = repmat(X_t, [1,1,xj_dim]);
Y_t        = repmat(Y_t, [1,1,yj_dim]);

% create distance grid and grid of hovering positions
D_mat_est  = ones(x_dim,y_dim, length(D_est));
Y_jt       = ones(x_dim,y_dim,yj_dim);
X_jt       = ones(x_dim,y_dim,xj_dim);
for k = 1:length(D_est)
    D_mat_est(:,:,k) = D_est(k);
    Y_jt(:,:,k)      = y_jt(k);
    X_jt(:,:,k)      = x_jt(k);
end

% calc the mle
P_xy = (X_t-X_jt).^2+(Y_t-Y_jt).^2 +params.sim.H.^2;
log_Pxy = log(P_xy);
numerator   = D_mat_est-sqrt(P_xy);
denominator = P_xy;
fraction = (numerator./denominator).^2;
% constants
P              = params.sim.P;
G_p            = params.sim.G_p;
beta_0         = params.sim.beta_0;
a              = params.sim.a;
sigma_0        = params.sim.sigma_0;
% final expression of the mle
factor   = (P*G_p*beta_0)./(a*sigma_0.^2);
func = sum(log_Pxy + factor*fraction,3);
% calc min 
[min_val,idx]   = min(func(:));
[idy,idx] = ind2sub(size(func),idx); 
end

