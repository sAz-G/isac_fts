%------------------------------------------------------------------------
% FUNCTION NAME: init_trajectory
% AUTHORS: Sharif Azem (sAz-G on GitHub) and Markus Krantzik (mardank on GitHub)
%
% DESCRIPTION: Calculates the initial trajectory. 
% STATUS:      Work in progress
% TODO:        Write code to handle situations, in which the initial
% trajectory lies out of the defined boundaries
% trajectory is out of the boundaries 
%
% INPUTS:
%   s_start - Starting point of the trajectory.
%   s_end   - Ending point of the trajectory.
%   N       - Number of trajectory points.
%   params  - Predefined parameters.
%
% OUTPUTS:
%   S_init  - Initial trajectory matrix.
%   V_init  - Initial velocity matrix.
%
% USAGE: [S_init, V_init] = init_trajectory(s_start, s_end, N, params)
%-----------------------------------------------------------------------

function [S_init, V_init] = init_trajectory(s_start, s_end, N, params)
    V_str = params.sim.V_str;
    T_f = params.sim.T_f;
    L_x = params.sim.L_x;
    L_y = params.sim.L_y;

    S_init = zeros(2, N);
    V_init = ones(2, N) .* V_str .* (s_end - s_start) ./ norm(s_end - s_start);

    S_init(:, 1:N) = s_start + V_str .* (1:N) .* (s_end - s_start) ./ norm(s_end - s_start) * T_f;
    
    %s_e = S_init(:, end); % Get last position of the initial trajectory
    
    %{
    if (s_e(1) > L_x) && (s_e(2) > L_y)  % Position out of bounds
        % Adjust the trajectory to stay within bounds
        if (abs(s_start(2) - L_y) > abs(s_start(1) - L_x))
            x_new = L_x;
            grad = s_start - s_end;
            grad = grad(2) ./ grad(1);
            y_new = s_e(2) + grad * (x_new - s_e(1));
            s_e(1) = x_new;
            s_e(2) = y_new;
            
        else
            y_new = L_y;
            grad_inv = s_start - s_end;
            grad_inv = grad_inv(1) ./ grad_inv(2);
            x_new = s_end(1) + grad_inv * (y_new - s_e(2));
            s_e(2) = y_new;
            s_e(1) = x_new;
        end

        % Recalculate the final trajectory
        V_str = norm(s_e - s_start) ./ ((N + 1) .* T_f);
        S_init(1, :) = linspace(s_start(1) + V_str * T_f, s_e(1), N);
        S_init(2, :) = linspace(s_start(2) + V_str * T_f, s_e(2), N);
    end 
    %}

end
