%------------------------------------------------------------------------
% FUNCTION NAME: init_trajectory
% AUTHOR: Sharif Azem     (TU-Darmstadt department 18, sAz-G on github)
%         Markus Krantzik (TU-Darmstadt department 18, mardank on github)
%
% DESCRIPTION: calculates the initial trajectory
%
% INPUTS:
%   s_start - starting point
%   s_end   - ending point
%   N       - amount of trajectory pionts
%   params - predefined parameters
%
% OUTPUTS:
%   S_init   - initial trajectory
%   V_init   - initial velocity 
%
% USAGE: [S_init, V_init] = init_trajectory(s_start, s_end,N, params)
%-----------------------------------------------------------------------

function [S_init, V_init] = init_trajectory(s_start, s_end,N, params)
V_str = params.sim.V_str;
T_f   = params.sim.T_f;
L_x = params.sim.L_x;
L_y = params.sim.L_y;

S_init = zeros(2,N);
V_init  = ones(2,N).*V_str.*(s_end - s_start)./norm(s_end - s_start);

S_init(:,1:N) = s_start + V_str.*(1:N).*(s_end - s_start)./norm(s_end - s_start)*T_f;

% check if the calculated trajectory is valid. if not, calculate a new one.
% (does not cover all cases but works fine)
s_e = S_init(:,end); % get last position at the initial trajectory

if (s_e(1) > L_x)&&(s_e(2) > L_y)                                    %  position out of bounds
    if (abs(s_start(2)-L_y) > abs(s_start(1)-L_x))
        x_new  = L_x;                                                    % calc new x
        s_e(1) = x_new;                             % calc new x
        grad   = s_start - s_end;                   % calc diff to obtain slope
        grad   = grad(2)./grad(1);                  % calc slope
        y_new  = s_e(2) + grad.*(x_new-s_e(1));   % calc y 
        s_e(2) = y_new;                             % set new y
    else
        y_new     = L_y;
        s_e(2)    = y_new;
        grad_inv  = s_start - s_end;
        grad_inv  = grad_inv(1)./grad_inv(2);
        x_new     = s_end(1) + grad_inv.*(y_new-s_e(2));
        s_e(1)    = x_new;
    end
    
    % calc final trajectory
    V_str = norm(s_e-s_start)./((N+1).*T_f);
    S_init(1,:) = linspace(s_start(1) + V_str.*T_f,s_e(1), N); % obtain linear trajectory. Note velocity is different than input
    S_init(2,:) = linspace(s_start(2) + V_str.*T_f,s_e(2), N);
elseif (s_e(1) > L_x)&&(s_e(2) < 0) % position out of bounds
    if (s_start(2) > abs(s_start(1)-L_x))
        x_new  = L_x;                                                    % calc new x
        s_e(1) = x_new;                             % calc new x
        grad   = s_start - s_end;                   % calc diff to obtain slope
        grad   = grad(2)./grad(1);                  % calc slop
        y_new  = s_e(2) + grad.*(x_new-s_e(1));   % calc y (no need to check boundary since x boundary is checked, and modification is linear)
        s_e(2) = y_new;                             % set new x
    else
        y_new     = 0;
        s_e(2)    = y_new;
        grad_inv  = s_start - s_end;
        grad_inv  = grad_inv(1)./grad_inv(2);
        x_new     = s_end(1) + grad_inv.*(y_new-s_e(2));
        s_e(1)    = x_new;
    end
    
    V_str = norm(s_e-s_start)./((N+1).*T_f);
    S_init(1,:) = linspace(s_start(1) + V_str.*T_f,s_e(1), N); % obtain linear trajectory. Note velocity is different than input
    S_init(2,:) = linspace(s_start(2) + V_str.*T_f,s_e(2), N);
elseif (s_e(1) < 0)&&(s_e(2) > L_y) % position out of bounds
    if (abs(s_start(2)-L_y) > s_start(1))
        x_new  = 0;                                                    % calc new x
        s_e(1) = x_new;                             % calc new x
        grad   = s_start - s_end;                   % calc diff to obtain slope
        grad   = grad(2)./grad(1);                  % calc slop
        y_new  = s_e(2) + grad.*(x_new-s_e(1));   % calc y (no need to check boundary since x boundary is checked, and modification is linear)
        s_e(2) = y_new;                             % set new x
    else
        y_new     = L_y;
        s_e(2)    = y_new;
        grad_inv  = s_start - s_end;
        grad_inv  = grad_inv(1)./grad_inv(2);
        x_new     = s_end(1) + grad_inv.*(y_new-s_e(2));
        s_e(1)    = x_new;
    end
    
    V_str = norm(s_e-s_start)./((N+1).*T_f);
    S_init(1,:) = linspace(s_start(1) + V_str.*T_f,s_e(1), N); % obtain linear trajectory. Note velocity is different than input
    S_init(2,:) = linspace(s_start(2) + V_str.*T_f,s_e(2), N);
elseif (s_e(1)  < 0)&&(s_e(2)  < 0) % position out of bounds
    if (s_start(2) > s_start(1))
        x_new  = 0;                
        % calc new x
        s_e(1) = x_new;                             % calc new x
        grad   = s_start - s_end;                   % calc diff to obtain slope
        grad   = grad(2)./grad(1);                  % calc slope
        y_new  = s_e(2) + grad.*(x_new-s_e(1));     % calc y 
        s_e(2) = y_new;                             % set new y
    else
        y_new     = 0;
        s_e(2)    = y_new;
        grad_inv  = s_start - s_end;
        grad_inv  = grad_inv(1)./grad_inv(2);
        x_new     = s_end(1) + grad_inv.*(y_new-s_e(2));
        s_e(1)    = x_new;
    end
    
    % final trajectory
    V_str = norm(s_e-s_start)./((N+1).*T_f);
    S_init(1,:) = linspace(s_start(1) + V_str.*T_f,s_e(1), N); % obtain linear trajectory. Note velocity is different than input
    S_init(2,:) = linspace(s_start(2) + V_str.*T_f,s_e(2), N);
elseif s_e(1) > L_x % x position out of bounds
    x_new  = L_x;                               % calc new x
    s_e(1) = x_new;                             % calc new x
    grad   = s_start - s_end;                   % calc diff to obtain slope
    grad   = grad(2)./grad(1);                  % calc slop
    y_new  = s_end(2) + grad.*(x_new-s_e(1));   % calc y 
    s_e(2) = y_new;                             % set new y
    
    % final trajectory
    V_str = norm(s_e-s_start)./((N+1).*T_f);
    S_init(1,:) = linspace(s_start(1) + V_str.*T_f,s_e(1), N);
    S_init(2,:) = linspace(s_start(2) + V_str.*T_f,s_e(2), N);
elseif s_e(1) < 0 % y position out of bounds
    x_new = 0;
    s_e(1) = x_new;
    grad   = s_start - s_end;
    grad   = grad(2)./grad(1);
    y_new  = s_end(2) + grad.*(x_new-s_e(1));
    s_e(2) = y_new;
    
    % final trajectory
    V_str = norm(s_e-s_start)./((N+1).*T_f); 
    S_init(1,:) = linspace(s_start(1) + V_str.*T_f,s_e(1), N);
    S_init(2,:) = linspace(s_start(2) + V_str.*T_f,s_e(2), N);
elseif s_e(2) > L_y  % y position out of bounds
    y_new     = L_y;
    s_e(2)    = y_new;
    grad_inv  = s_start - s_end;
    grad_inv  = grad_inv(1)./grad_inv(2);
    x_new     = s_end(1) + grad_inv.*(y_new-s_e(2));
    s_e(1)    = x_new;
    
    % final trajectory
    V_str = norm(s_e-s_start)./((N+1).*T_f);
    S_init(1,:) = linspace(s_start(1) + V_str.*T_f,s_e(1), N);
    S_init(2,:) = linspace(s_start(2) + V_str.*T_f,s_e(2), N);
elseif s_e(2) < 0  % y position out of bounds
    y_new    = 0;
    s_e(2)   = y_new;
    grad_inv = s_start - s_end;
    grad_inv = grad_inv(1)./grad_inv(2);
    x_new    = s_end(1) + grad_inv.*(y_new-s_e(2));
    s_e(1)   = x_new;
    
    %final trajectory
    V_str = norm(s_e-s_start)./((N+1).*T_f);
    S_init(1,:) = linspace(s_start(1) + V_str.*T_f,s_e(1), N);
    S_init(2,:) = linspace(s_start(2) + V_str.*T_f,s_e(2), N);
end

end

