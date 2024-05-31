%------------------------------------------------------------------------
% FUNCTION NAME: optimize_m
% AUTHOR: Sharif Azem (sAz-G on github), Markus Krantzik (mardank on github)
%
%
% DESCRIPTION:  optimize_m optimizes the trajectory for given initial values at the
% mth stage. This function returns the final results. 
%
% INPUTS:
% 
% E_m               - the available energy at the mth stage. This is needed for the constraint. 
% s_c               - position of the communication user 
% S_hover           - all hover points including the mth stage. The hover points K
%                     hover points of the mth stage are going to be optimized.
% S_total           - all trajectory points including these of the mth stage. The
%                     points of the mth stage are going to be optimized.
% params            - constants that are defined before running the whole program
%
% s_s               - starting point
% s_b               - base station position
% s_t               - (estimated) target position 
%
% OUTPUTS:
%   S_opt     - optimal trajectory of stage m
%   V_opt     - Velocity matrix
%   xi_opt    - opt. variable xi
%   dlta_opt - opt. variable dlta
%   CRB_opt   - opt. crb used in the objective function J
%   R_opt     - opt. rate function used in the objective function J
%   CRB_iter  - value of CRB_opt at each iteration
%   R_iter    - value of R_opt at each iteration  
%   J         - objective function 
%
% USAGE: [S_opt, V_opt, xi_opt, dlta_opt, CRB_opt, R_opt, CRB_iter, R_iter , J] = optimize_m(E_m, s_c, S_hover, S_total,s_t,s_s, params)
%
%-----------------------------------------------------------------------

function [S_opt, V_opt, xi_opt, dlta_opt, CRB_opt, R_opt, CRB_iter, R_iter , J] = optimize_m(E_m, s_c, S_hover, S_total,s_t,s_s, params)

% constants 
L_x   = params.sim.L_x;     % boundary in the x direction 
L_y   = params.sim.L_y;     % boundary in the y direction 
V_max = params.sim.V_max;   % max possibe velocity

v_0     = params.energy.v_0;    % energy parameter 
 
iter    = params.sim.iter;    % amount of iterations for the optimization algo
mu      = params.sim.mu;      % sense target every mu points. 1,2,mu,...,N_m
eta     = params.sim.eta;     % trade of constant for the object function 
w_star  = params.sim.w_star;  % step size
N_stg   = params.sim.N_stg;   % amount of trajectory points at each stage 
K_stg   = params.sim.K_stg;   % amount of trajectory hover points at each stage 
best_val = 1e13;

%% initial values 
S_init = S_total(:,end-N_stg+1:end); 
V_init = nan(size(S_init));
V_init(:,1)       = (S_init(:,1)-s_s)./params.sim.T_f;
V_init(:,2:end)   = (S_init(:,2:end)-S_init(:,1:end-1))./params.sim.T_f;
dlta_square_init = sqrt(1 + norms(V_init, 2, 1).^4/(4*params.energy.v_0^4))...
                        - norms(V_init, 2 ,1).^2/(2*params.energy.v_0^2);
S_init_valid = S_init;
S_valid      = S_init*1.1; % initial value for the valid solution. used only in case of nan in the first iteration  
N_nan        = 0;             % amount of consecutive nan solutions

J = nan(iter,1);
CRB_iter =  nan(iter,1);
R_iter = nan(iter,1);

for u = 1:iter  % optimization loop

% assign optimized hover points and then optimized trajectory
S_hover(:,end-K_stg+1:end)  = [S_init(1, mu:mu:N_stg); S_init(2, mu:mu:N_stg)]; 
S_total(:,end-N_stg+1:end)  = S_init;

% gradient of crb in x direction 
crb_grad_x = crb_grad(S_hover, s_t, params, 'x');

% gradient of crb in y direction 
crb_grad_y = crb_grad(S_hover, s_t, params, 'y'); 
 
% gradient of rate in x direction 
rate_grad_x = rate_grad(S_total(:,end-N_stg+1:end), s_c, params, 'x',size(S_total,2));

% gradient of rate in y direction 
rate_grad_y = rate_grad(S_total(:,end-N_stg+1:end), s_c, params, 'y',size(S_total,2));

cvx_begin % cvx entry point
    if strcmp(params.opt_settings.solver, 'mosek')
        cvx_solver mosek         % use solver mosek
    elseif strcmp(params.opt_settings.solver, 'sedumi')
            cvx_solver sedumi    % use solver sedumi
    elseif strcmp(params.opt_settings.solver, 'sdpt3')
            cvx_solver sdpt3     % use solver sdpt3
    end

    % set the precision of the chosen solver
    if strcmp(params.opt_settings.precision, 'low')
            cvx_precision low  
    elseif strcmp(params.opt_settings.precision, 'medium')
            cvx_precision medium  
    elseif strcmp(params.opt_settings.precision, 'high')
            cvx_precision high  
    elseif strcmp(params.opt_settings.precision, 'best')
            cvx_precision best
    else
            cvx_precision default  
    end

    if ~params.opt_settings.screen_out
        cvx_quiet true
    end
    
    % optimization variables
    variable S(2,N_stg)     % trajectory 
    variable V(2,N_stg)     % velocity
    variable dlta(1,N_stg) % dlta (see paper)
    variable xi(1,N_stg)    % xi (see paper)
    
    % first degree taylor series expansion of crb in x direction 
    crb_taylor_x = sum(crb_grad_x.*(S(1,mu:mu:N_stg) - S_init(1, mu:mu:N_stg)));
    % first degree taylor series expansion of crb in y direction
    crb_taylor_y = sum(crb_grad_y.*(S(2,mu:mu:N_stg) - S_init(2, mu:mu:N_stg))); 
    
    % first degree taylor series expansion of crb 
    CRB_taylor =  crb_taylor_x + crb_taylor_y;

    % first degree of taylor series expansion in x direction
    rate_taylor_x = sum(rate_grad_x.*(S(1,:) - S_init(1,:))); 
    % first degree of taylor series expansion in y direction
    rate_taylor_y = sum(rate_grad_y.*(S(2,:) -  S_init(2,:))); 
    
    % first degree taylor series expansion of the rate
    R_taylor =  rate_taylor_x + rate_taylor_y;

    % minimization
    minimize(eta.*(CRB_taylor)-(1-eta).*R_taylor./(params.sim.B));
    
    % constraints    
    subject to
        E_m     >= calc_constraint_energy(V,dlta, params); % energy constraint
        (S(:,1) - s_s)./params.sim.T_f == V(:,1); % definition of the velocity
        (S(:,2:end) - S(:,1:end-1))./params.sim.T_f == V(:,2:end); % definition of the velocity
        V_max   >= norms(V, 2, 1);    % velocity constraint
        dlta   >= 0;                 % dlta is not negative 
        xi      >= 0;                 % xi is not negative
        S(1,:)  >= 0;                 % all x points are not negative                                                      
        S(2,:)  >= 0;                 % all y points are not negative
        L_x     >= S(1,:);            % all x points are bounded by L_x
        L_y     >= S(2,:);            % all y points are bounded by L_y
        
        for i = 1:N_stg % add constraints 51a and 51b from the paper
           (norm(V_init(:, i))./v_0).^2 + (2./v_0.^2).*V_init(:, i).'*(V(:,i)-V_init(:,i)) >= pow_pos(inv_pos(dlta(i)), 2)-xi(i);
           dlta_square_init(i) + 2.*sqrt(dlta_square_init(i)).*(dlta(i) - sqrt(dlta_square_init(i))) >=  xi(i);
        end
cvx_end % end of cvx calculations 

% handle solution
if isnan(cvx_optval) %|| isinf(cvx_optval)% handle nan solutions, set new starting point and override invalid solutions with old ones
    S      = S_valid;            % use the last valid solution, override current solution 
    S_init = S_init_valid;       % use the last S_init which led to a valid solution 
    N_nan  = N_nan +1;
    if u < iter % set an initial trajectory for the next iteration  
        %S_init            = S_init + ((1-w_star).^(N_nan)).*(S-S_init); % choose a point far away from the last initial point
        S_init            = S_init + (0.97.^(N_nan.^N_nan)).*(S-S_init); % start search from 0.97 to 0 to find a point that leads to a valid solution
        %S_init            = S_init + w_star.*(S-S_init);

        V_init(:,1)       = (S_init(:,1)-s_s)./params.sim.T_f;
        V_init(:,2:end)   = (S_init(:,2:end)-S_init(:,1:end-1))./params.sim.T_f;
        dlta_square_init = sqrt(1 + norms(V_init, 2, 1).^4/(4*params.energy.v_0^4))...
                                - norms(V_init, 2 ,1).^2/(2*params.energy.v_0^2);
    else % if nan solution and last iteration, assign last valid solutions
        break; % stop optimizing if reached last iteration
    end
else% if the solution is valid (not nan)
    N_nan = 0; % set this to 0 everytime the program provides a valid solution
    %store valid optimal solution in buffers. These are used only if in the
    %next iteration the algo provides nan solution 
    if abs(best_val) > abs(cvx_optval)
        S_valid         = S;
        S_init_valid    = S_init; 
        V_valid         = V;
        dlta_opt_valid = dlta;
        xi_valid        = xi;
        R_opt_valid     = R_taylor;
        CRB_opt_valid   = CRB_taylor;
        best_val        = cvx_optval;        
    end
   
    J(u) = eta.*(CRB_taylor)-(1-eta).*R_taylor./(params.sim.B);
    R_iter(u)    = R_taylor;
    CRB_iter(u)  = CRB_taylor;
    
    if (u < iter) && (~( abs(cvx_optval) < params.opt_settings.threshold)) % set an initial values for the next iteration  
        S_init            = S_init + w_star.*(S-S_init);
        V_init(:,1)       = (S_init(:,1)-s_s)./params.sim.T_f;
        V_init(:,2:end)   = (S_init(:,2:end)-S_init(:,1:end-1))./params.sim.T_f;
        dlta_square_init = sqrt(1 + norms(V_init, 2, 1).^4/(4*params.energy.v_0^4))...
                                - norms(V_init, 2 ,1).^2/(2*params.energy.v_0^2);
    else
        break; % stop optimizing
    end


end % end of if isnan(cvx_optval)
end % end of for loop and optimization

S_opt       = S_valid;              % store trajectory results of each iteration 
V_opt       = V_valid;              % store velocity results of each iteration 
dlta_opt   = dlta_opt_valid;      % store detla results of each iteration 
xi_opt      = xi_valid;             % store xi results of each iteration 
R_opt       = R_opt_valid;          % store rate results of each iteration 
CRB_opt     = CRB_opt_valid;        % store crb results  of each iteration 

end

