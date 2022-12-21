function [idx, idy] = get_min(D_est,x_t,y_t,x_jt,y_jt,H)
% get_min(est, params) is the function, for which we want to find the
% minimum(objective function (of)).
% 
% Note that this is not the original of, which is the log of a gaussian
% distribution (eq. 58 in the paper). The gaussian distribution is defined
% as eq. 57 in the paper.
%  
% After reshaping the log of eq. 57, and then modifing eq. 58, we get the
% following expression:
%
% min(L(f)) = 
% sum(log((x_{j}^{m}-x_{t})^2+(y_{j}^{m}-y_{t})^2))
% + sum(frac{dhat_{s}^{m}-sqrt((x_{j}^{m}-x_{t})^2+(y_{j}^{m}-y_{t})^2))}{(x_{j}^{m}-x_{t})^2+(y_{j}^{m}-y_{t})^2)})
%
% Arguments:
% x_t,y_t: With respect to this parameters we minimum the of. These are the
% possible target coordinates.
% D_est: The vector of estimated positions up to the current estimation.

x_dim = length(x_t); % dimension in x direction of the matrices after using meshgrid
y_dim = length(y_t); % dimension in y direction of the matrices after using meshgrid
xj_dim = length(y_jt); % dimension in x direction of the matrices after using meshgrid
yj_dim = length(x_jt); % dimension in y direction of the matrices after using meshgrid

[X_t, Y_t] = meshgrid(y_t, x_t);
X_t        = repmat(X_t, [1,1,xj_dim]);
Y_t        = repmat(Y_t, [1,1,yj_dim]);

D_mat_est  = ones(x_dim,y_dim, length(D_est));
Y_jt       = ones(x_dim,y_dim,yj_dim);
X_jt       = ones(x_dim,y_dim,xj_dim);
H_mat      = ones(x_dim,y_dim, length(D_est)).*H.^2;

for k = 1:length(D_est)
    D_mat_est(:,:,k) = D_est(k);
    Y_jt(:,:,k)      = x_jt(k);
    X_jt(:,:,k)      = y_jt(k);
end

P_xy = (X_t-X_jt).^2+(Y_t-Y_jt).^2 + H_mat.^2;
log_Pxy = log(P_xy);

numerator   = D_mat_est-sqrt(P_xy);
denominator = P_xy;

fraction = (numerator./denominator).^2;

func = sum(log_Pxy + fraction,3);

[~,idx]   = min(func(:));
[idx,idy] = ind2sub(size(func),idx); 


end

