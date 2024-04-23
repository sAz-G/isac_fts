%------------------------------------------------------------------------
% FUNCTION NAME: optimize_m
% AUTHORS: Sharif Azem (sAz-G on GitHub) and Markus Krantzik (mardank on GitHub)
%
% DESCRIPTION: Optimizes the trajectory for given initial values at the mth stage.
%
% INPUTS:
%   E_m         - Available energy at the mth stage (needed for the constraint).
%   s_c         - Position of the communication user.
%   S_hover     - All hover points including the mth stage. The hover points K
%                 hover points of the mth stage are going to be optimized.
%   S_total     - All trajectory points including those of the mth stage. The
%                 points of the mth stage are going to be optimized.
%   s_t         - (Estimated) target position.
%   s_s         - Starting point.
%   params      - Constants defined before running the whole program.
%
% OUTPUTS:
%   S_opt       - Optimal trajectory of stage m.
%   V_opt       - Velocity matrix.
%   xi_opt      - Optimal variable xi.
%   delta_opt   - Optimal variable delta.
%   CRB_opt     - Optimal CRB used in the objective function J.
%   R_opt       - Optimal rate function used in the objective function J.
%   CRB_iter    - Value of CRB_opt at each iteration.
%   R_iter      - Value of R_opt at each iteration.
%   J           - Objective function.
%
% USAGE: [S_opt, V_opt, xi_opt, delta_opt, CRB_opt, R_opt, CRB_iter, R_iter, J] = optimize_m(E_m, s_c, S_hover, S_total, s_t, s_s, params)
%
%-----------------------------------------------------------------------

function [S_opt, V_opt, xi_opt, delta_opt, CRB_opt, R_opt, CRB_iter, R_iter, J] = optimize_m(E_m, s_c, S_hover, S_total, s_t, s_s, params)

% Constants 
L_x = params.sim.L_x;       % Boundary in the x direction.
L_y = params.sim.L_y;       % Boundary in the y direction.
V_max = params.sim.V_max;   % Maximum possible velocity.
v_0 = params.energy.v_0;    % Energy parameter.
iter = params.sim.iter;     % Number of iterations for the optimization algorithm.
mu = params.sim.mu;         % Sense target every mu points: 1, 2, mu, ..., N_m.
eta = params.sim.eta;       % Trade-off constant for the objective function.
w_star = params.sim.w_star; % Step size.
N_stg = params.sim.N_stg;   % Number of trajectory points at each stage.
K_stg = params.sim.K_stg;   % Number of trajectory hover points at each stage.
best_val = 1e13;            % Initial best value.

%% Initial values 
S_init = S_total(:, end - N_stg + 1:end);
V_init = nan(size(S_init));
V_init(:,1) = (S_init(:,1) - s_s) ./ params.sim.T_f;
V_init(:,2:end) = (S_init(:,2:end) - S_init(:,1:end-1)) ./ params.sim.T_f;
delta_square_init = sqrt(1 + norms(V_init, 2, 1).^4 / (4 * params.energy.v_0^4)) - norms(V_init, 2 ,1).^2 / (2 * params.energy.v_0^2);
S_init_valid = S_init;
S_valid = S_init * 1.1; % Initial value for the valid solution. Used only in case of NaN in the first iteration.
N_nan = 0; % Counter for consecutive NaN solutions.

J = nan(iter, 1);
CRB_iter = nan(iter, 1);
R_iter = nan(iter, 1);

for u = 1:iter % Optimization loop

    % Assign optimized hover points and then optimized trajectory
    S_hover(:, end - K_stg + 1:end) = [S_init(1, mu:mu:N_stg); S_init(2, mu:mu:N_stg)];
    S_total(:, end - N_stg + 1:end) = S_init;

    % Gradient of CRB in x direction
    crb_grad_x = crb_grad(S_hover, s_t, params, 'x');

    % Gradient of CRB in y direction
    crb_grad_y = crb_grad(S_hover, s_t, params, 'y');

    % Gradient of rate in x direction
    rate_grad_x = rate_grad(S_total(:, end - N_stg + 1:end), s_c, params, 'x', size(S_total, 2));

    % Gradient of rate in y direction
    rate_grad_y = rate_grad(S_total(:, end - N_stg + 1:end), s_c, params, 'y', size(S_total, 2));

    cvx_begin % CVX entry point
        if strcmp(params.opt_settings.solver, 'mosek')
            cvx_solver mosek % Use solver Mosek.
        elseif strcmp(params.opt_settings.solver, 'sedumi')
            cvx_solver sedumi % Use solver Sedumi.
        elseif strcmp(params.opt_settings.solver, 'sdpt3')
            cvx_solver sdpt3 % Use solver SDPT3.
        end

        % Set the precision of the chosen solver
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

        % Optimization variables
        variable S(2, N_stg) % Trajectory
        variable V(2, N_stg) % Velocity
        variable delta(1, N_stg) % Delta (see paper)
        variable xi(1, N_stg) % Xi (see paper)

        % First degree Taylor series expansion of CRB in x direction
        crb_taylor_x = sum(crb_grad_x .* (S(1, mu:mu:N_stg) - S_init(1, mu:mu:N_stg)));

        % First degree Taylor series expansion of CRB in y direction
        crb_taylor_y = sum(crb_grad_y .* (S(2, mu:mu:N_stg) - S_init(2, mu:mu:N_stg)));

        % First degree Taylor series expansion of CRB
        CRB_taylor = crb_taylor_x + crb_taylor_y;

        % First degree of Taylor series expansion in x direction
        rate_taylor_x = sum(rate_grad_x .* (S(1, :) - S_init(1, :)));

        % First degree of Taylor series expansion in y direction
        rate_taylor_y = sum(rate_grad_y .* (S(2, :) - S_init(2, :)));

        % First degree Taylor series expansion of the rate
        R_taylor = rate_taylor_x + rate_taylor_y;

        % Minimization
        minimize(eta * CRB_taylor - (1 - eta) * R_taylor / params.sim.B);

        % Constraints    
        subject to
            E_m >= calc_constraint_energy(V, delta, params); % Energy constraint
            (S(:, 1) - s_s) ./ params.sim.T_f == V(:, 1); % Definition of the velocity
            (S(:, 2:end) - S(:, 1:end-1)) ./ params.sim.T_f == V(:, 2:end); % Definition of the velocity
            V_max >= norms(V, 2, 1); % Velocity constraint
            delta >= 0; % Delta is not negative
            xi >= 0; % Xi is not negative
            S(1, :) >= 0; % All x points are non-negative                                                      
            S(2, :) >= 0; % All y points are non-negative
            L_x >= S(1, :); % All x points are bounded by L_x
            L_y >= S(2, :); % All y points are bounded by L_y
            
            % Add constraints 51a and 51b from the paper
            for i = 1:N_stg
                (norm(V_init(:, i)) / v_0)^2 + (2 / v_0^2) * V_init(:, i).' * (V(:, i) - V_init(:, i)) >= pow_pos(inv_pos(delta(i)), 2) - xi(i);
                delta_square_init(i) + 2 * sqrt(delta_square_init(i)) * (delta(i) - sqrt(delta_square_init(i))) >=  xi(i);
            end
    cvx_end % End of CVX calculations

    % Handle solution
    if isnan(cvx_optval) % Handle NaN solutions
        S = S_valid; % Use the last valid solution, override current solution
        S_init = S_init_valid; % Use the last S_init which led to a valid solution 
        N_nan = N_nan + 1;
        if u < iter % Set initial trajectory for the next iteration  
            S_init = S_init + (0.97^(N_nan^N_nan)) * (S - S_init); % Start search from 0.97 to 0 to find a point that leads to a valid solution
            V_init(:, 1) = (S_init(:, 1) - s_s) ./ params.sim.T_f;
            V_init(:, 2:end) = (S_init(:, 2:end) - S_init(:, 1:end-1)) ./ params.sim.T_f;
            delta_square_init = sqrt(1 + norms(V_init, 2, 1).^4 / (4 * params.energy.v_0^4)) - norms(V_init, 2 ,1).^2 / (2 * params.energy.v_0^2);
        else % If NaN solution and last iteration, assign last valid solutions
            break; % Stop optimizing if reached last iteration
        end
    else % If the solution is valid (not NaN)
        N_nan = 0; % Set this to 0 every time the program provides a valid solution
        % Store valid optimal solution in buffers
        if abs(best_val) > abs(cvx_optval)
            S_valid = S;
            S_init_valid = S_init; 
            V_valid = V;
            delta_opt_valid = delta;
            xi_valid = xi;
            R_opt_valid = R_taylor;
            CRB_opt_valid = CRB_taylor;
            best_val = cvx_optval;        
        end

        J(u) = eta * CRB_taylor - (1 - eta) * R_taylor / params.sim.B;
        R_iter(u) = R_taylor;
        CRB_iter(u) = CRB_taylor;

        if (u < iter) && (~ (abs(cvx_optval) < params.opt_settings.threshold)) % Set initial values for the next iteration  
            S_init = S_init + w_star * (S - S_init);
            V_init(:, 1) = (S_init(:, 1) - s_s) ./ params.sim.T_f;
            V_init(:, 2:end) = (S_init(:, 2:end) - S_init(:, 1:end-1)) ./ params.sim.T_f;
            delta_square_init = sqrt(1 + norms(V_init, 2, 1).^4 / (4 * params.energy.v_0^4)) - norms(V_init, 2 ,1).^2 / (2 * params.energy.v_0^2);
        else
            break; % Stop optimizing
        end
    end % End of if isnan(cvx_optval)
end % End of for loop and optimization

% Store results
S_opt = S_valid; % Optimal trajectory results of each iteration 
V_opt = V_valid; % Optimal velocity results of each iteration 
delta_opt = delta_opt_valid; % Optimal delta results of each iteration 
xi_opt = xi_valid; % Optimal xi results of each iteration 
R_opt = R_opt_valid; % Optimal rate results of each iteration 
CRB_opt = CRB_opt_valid; % Optimal CRB results of each iteration 

end
