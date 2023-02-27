function [S_m, V_m, xi_m, delta_m, CRB_iter, R_iter] = optimize_m_debug(E_m, s_c,S_hover, S_total, delta_square_init,s_t,s_s, V_init, params)
% optimize_m_debug optimizes the trajectory for given initial values at the
% mth stage
%   
% input: 
% 
% E_m               - the available energy at the mth stage. This is needed for the constraint. 
% s_c               - position of the communication user 
% S_hover           - all hover points including the mth stage. The hover points K
%                     hover points of the mth stage are going to be optimized.
% S_total           - all trajectory points including these of the mth stage. The
%                     points of the mth stage are going to be optimized.
% delta_square_init - initial delta square value passed to the function. It changes in every iteration
%                     of the optimization algorithm 
% V_init            - initial velocity 
% params            - constants that are defined before running the whole program
%
% return:
% optimization variables of each iteration 

% constants 
L_x   = params.sim.L_x;     % boundary in the x direction 
L_y   = params.sim.L_y;     % boundary in the y direction 
V_max = params.sim.V_max;   % max possibe velocity

v_0     = params.energy.v_0;    % energy parameter 

iter    = params.sim.iter;   % amount of iterations for the optimization algo
mu      = params.sim.mu;     % sense target every mu points. 1,2,mu,...,N_m
eta     = params.sim.eta;    % trade of constant for the object function 
w_star  = params.sim.w_star; % step size
N_stg   = params.sim.N_stg;  % amount of trajectory points at each stage 
K_stg   = params.sim.K_stg;   % amount of trajectory hover points at each stage 

CRB_iter = nan(1,iter);       % store crb values of each iteration 
R_iter   = nan(1,iter);       % store rate values of each iteration 

%% Debug parameter 

S_mat     = zeros(2, N_stg, iter); % store trajectory results of each iteration
V_mat     = zeros(2, N_stg, iter); % store velocity results of each iteration
delta_mat = zeros(N_stg, iter);    % store delta results of each iteration 
xi_mat    = zeros(N_stg, iter);    % store xi results of each iteration 

S_init = S_total(:,end-N_stg+1:end); % initial trajectory points. These points are going to be optimized

for u = 1:iter  % optimization loop

% assign optimized hover points and then optimized trajectory
S_hover(:,end-K_stg+1:end)  = [S_init(1, mu:mu:N_stg); S_init(2, mu:mu:N_stg)]; 
S_total(:,end-N_stg+1:end)  = S_init;

% gradient of crb in x direction 
crb_grad_x = crb_grad(S_hover, s_t, params, 'x');

% gradient of crb in y direction 
crb_grad_y = crb_grad(S_hover, s_t, params, 'y'); 
 
% gradient of rate in x direction 
rate_grad_x = rate_grad(S_total, s_c, params, 'x');

% gradient of rate in y direction 
rate_grad_y = rate_grad(S_total, s_c, params, 'y');

cvx_begin % cvx entry point
    cvx_solver mosek    % use solver mosek
    cvx_precision high  
    
    % optimization variables
    variable S(2,N_stg)     % trajectory 
    variable V(2,N_stg)     % velocity
    variable delta(1,N_stg) % delta (see paper)
    variable xi(1,N_stg)    % xi (see paper)
    
    % first degree taylor series expansion of crb in x direction 
    crb_taylor_x = sum(crb_grad_x.*(S(1,mu:mu:N_stg) - S_init(1, mu:mu:N_stg)));
    % first degree taylor series expansion of crb in y direction
    crb_taylor_y = sum(crb_grad_y.*(S(2,mu:mu:N_stg) - S_init(2, mu:mu:N_stg))); 
    
    % first degree taylor series expansion of crb 
    CRB_taylor =  crb_taylor_x + crb_taylor_y;

    % first degree of taylor series expansion in x direction
    rate_taylor_x = sum(rate_grad_x .* (S(1,:) - S_init(1,:))); 
    % first degree of taylor series expansion in y direction
    rate_taylor_y = sum(rate_grad_y .* (S(2,:) -  S_init(2,:))); 
    
    % first degree taylor series expansion of the rate
    R_taylor =  rate_taylor_x + rate_taylor_y;
 

    %% Objective function
    % minimization
    minimize(eta.*CRB_taylor-(1-eta).*R_taylor);

    %% constraints
    
    subject to
        E_m     >= calc_constraint_energy(V,delta, params); % energy constraint
        V_max   >= norms(V, 2, 1);    % velocity constraint
        (S(:,1) - s_s)./params.sim.T_f == V(:,1);
        (S(:,2:end) - S(:,1:end-1))./params.sim.T_f == V(:,2:end);
        delta   >= 0;                 % delta is not negative 
        xi      >= 0;                 % xi is not negative
        S(1,:)  >= 0;                 % all x points are not negative                                                      
        S(2,:)  >= 0;                 % all y points are not negative
        L_x     >= S(1,:);            % all x points are bounded by L_x
        L_y     >= S(2,:);            % all y points are bounded by L_y
       
        for i = 1:N_stg % add constraints 51a and 51b from the paper
           (norm(V_init(:, i))./v_0).^2 + (2./v_0.^2).*V_init(:, i).'*(V(:,i)-V_init(:,i)) >= pow_pos(inv_pos(delta(i)), 2)-xi(i);
           delta_square_init(i) + 2.*sqrt(delta_square_init(i)).*(delta(i) - sqrt(delta_square_init(i))) >=  xi(i);
        end
        
cvx_end % end of cvx calculations 

if u < iter % set an initial trajectory for the next iteration  
    S_init            = S_init + w_star.*(S-S_init);
    V_init(:,1)       = (S_init(:,1)-s_s)./params.sim.T_f;
    V_init(:,2:end)   = (S_init(:,2:end)-S_init(:,1:end-1))./params.sim.T_f;
    delta_square_init = sqrt(1 + norms(V_init, 2, 1).^4/(4*params.energy.v_0^4))...
                            - norms(V_init, 2 ,1).^2/(2*params.energy.v_0^2);
else % set S_init to the cvx result at the last iteration 
    S_init = S;
end

%V_init = V;                             % set initial velocity for the next iteration 
%delta_square_init = pow_pos(delta,2);   % set initial delta_square for the next iteration

S_mat(:, :, u) = S_init;      % store trajectory results of each iteration 
V_mat(:, :, u) = V;           % store velocity results of each iteration 
delta_mat(:, u) = delta;   % store detla results of each iteration 
xi_mat(:, u) = xi;         % store xi results of each iteration 

R_iter(u)    = R_taylor;         % store rate results of each iteration 
CRB_iter(u)  = CRB_taylor;       % store crb results  of each iteration 
end % end of for loop and optimization

% return output variables. This is the final optimization result
S_m         = S_mat;      % optimal trajectory 
V_m         = V_mat;      % optimal velocity
xi_m        = xi_mat;     % optimal xi
delta_m     = delta_mat;  % optimal delta
end

